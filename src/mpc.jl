function init_x(m, circuit, demand; opt::Option=Option())
    num_buses = length(circuit.bus)
    num_gens = length(circuit.gen)
    T = size(demand.pd,2)

    if !isempty(opt.neg_g)
        active_gen = findall(x -> !(x in opt.neg_g), collect(1:num_gens))
    else
        active_gen = collect(1:num_gens)
    end

    Vm = zeros(T, num_buses)
    for b=1:num_buses
        Vm[:,b] .= 0.5*(circuit.bus[b].Vmax + circuit.bus[b].Vmin)
    end
    Va = circuit.bus[circuit.busref].Va * ones(T, num_buses)

    setvalue(m[:Vm], Vm)
    setvalue(m[:Va], Va)

    Pg = zeros(T, num_gens)
    Qg = zeros(T, num_gens)

    for g=1:num_gens
        Pg[:,g] .= 0.5*(circuit.gen[g].Pmax + circuit.gen[g].Pmin)
        Qg[:,g] .= 0.5*(circuit.gen[g].Qmax + circuit.gen[g].Qmin)
    end

    setvalue(m[:Pg], Pg)
    setvalue(m[:Qg], Qg)

    if opt.sc_constr
        for g=1:num_gens
            setvalue(m[:Pg_f][g], Pg[1,g])
            setvalue(m[:Qg_f][g], Qg[1,g])
        end

        for b=1:num_buses
            setvalue(m[:Vm_f][b], Vm[1,b])
            setvalue(m[:Va_f][b], Va[1,b])
        end
    end

    if opt.freq_ctrl
        for t=1:T, g in active_gen
            setvalue(m[:Pg_ref][t,g], Pg[g])
        end
        setvalue(m[:omega], zeros(T))
    end

    if opt.load_shed
        setvalue(m[:sigma_p], zeros(T, num_buses))
        setvalue(m[:sigma_m], zeros(T, num_buses))
    end
end

function get_profname(case, T, load_scale, ramp_scale, warm)
    profname = "ipopt_"*splitext(basename(case))[1]
    profname = profname*"_t"*string(T)
    profname = profname*"_ls"*string(load_scale)*"_rs"*string(ramp_scale)

    if warm == "cold"
    profname = profname*"_cold_phase0"
    elseif warm == "shift_copy"
    profname = profname*"_warm_phase0"
    else
    profname = profname*"_warm_phase1"
    end

    return profname
end

function get_largest_gen(m_cur, circuit; howmany=1)
    Pg = getvalue(m_cur[:Pg])[1,:]
    idx = sortperm(Pg, rev=true)
    total = 0

    println("Top-", howmany, " generator:")
    for i=1:howmany
        println("Gen id: ", idx[i], " value: ", Pg[idx[i]])
        total += Pg[idx[i]]
    end
    println("Total: ", total)

    return idx[1:howmany]
end

function copyvars(inner_to, to_start, inner_from, from_start, n)
    copyto!(inner_to.x, to_start, inner_from.x, from_start, n)
    copyto!(inner_to.mult_x_L, to_start, inner_from.mult_x_L, from_start, n)
    copyto!(inner_to.mult_x_U, to_start, inner_from.mult_x_U, from_start, n)
end

function copyconstrs(inner_to, to_start, inner_from, from_start, n)
    copyto!(inner_to.g, to_start, inner_from.g, from_start, n)
    copyto!(inner_to.mult_g, to_start, inner_from.mult_g, from_start, n)
end

function get_cutgen(m, circuit, T)
    Pg = getvalue(m[:Pg])
    active = []

    if length(circuit.bus) in [73, 300]
        active = findall(x -> x > 0 && x <= 1, Pg[1,:])
    else
        active = findall(x -> x > 0, Pg[1,:])
    end

    order = sortperm(Pg[1,active], rev=true)
    max_g = argmax(Pg[1,active])
    min_g = argmin(Pg[1,active])
    med_g = order[max(1, floor(Int, 0.02*length(order)))]

    @printf("Generataor (max) %d: %.2f MW\n", active[max_g], Pg[1,active[max_g]]*circuit.baseMVA)
    @printf("Generataor (min) %d: %.2f MW\n", active[min_g], Pg[1,active[min_g]]*circuit.baseMVA)
    @printf("Generataor (med) %d: %.2f MW\n", active[med_g], Pg[1,active[med_g]]*circuit.baseMVA)

   return active[max_g]
end

function get_cutline(m, circuit, T)
    num_buses = length(circuit.bus)

    line = circuit.line
    ybus = circuit.ybus
    yline = circuit.yline
    busdict = circuit.busdict
    frombus = circuit.frombus
    tobus = circuit.tobus

    fromflow = zeros(T, length(line))
    toflow = zeros(T, length(line))

    Vm = getvalue(m[:Vm])
    Va = getvalue(m[:Va])

    for t=1:T,b=1:num_buses
        for l in frombus[b]
            fromflow[t,l] = (yline[l].YffR + ybus[b].YshR)*Vm[t,b]^2 +
                Vm[t,b]*Vm[t,busdict[line[l].to]]*
            (yline[l].YftR*cos(Va[t,b]-Va[t,busdict[line[l].to]]) +
             yline[l].YftI*sin(Va[t,b]-Va[t,busdict[line[l].to]]))
        end

        for l in tobus[b]
            toflow[t,l] = (yline[l].YttR + ybus[b].YshR)*Vm[t,b]^2 +
                Vm[t,b]*Vm[t,busdict[line[l].from]]*
            (yline[l].YtfR*cos(Va[t,b]-Va[t,busdict[line[l].from]]) +
             yline[l].YftI*sin(Va[t,b]-Va[t,busdict[line[l].from]]))
        end
    end

    abs_fromflow = abs.(fromflow)
    idx1 = findall(abs_fromflow .> 0)
    r1 = argmin(abs_fromflow[idx1])
    @printf("Line (%5d,%5d) has absolute minimum from-flow %10.8f\n",
            line[idx1[r1][2]].from, line[idx1[r1][2]].to,
            abs_fromflow[idx1[r1]])

    abs_toflow = abs.(toflow)
    idx2 = findall(abs_toflow .> 0)
    r2 = argmin(abs_toflow[idx2])
    @printf("Line (%5d,%5d) has absolute minimum to-flow %10.8f\n",
            line[idx2[r2][2]].from, line[idx2[r2][2]].to,
            abs_toflow[idx2[r2]])

    if abs_fromflow[idx1[r1]] < abs_toflow[idx2[r2]]
        return idx1[r1][2]
    else
        return idx2[r2][2]
    end
end

function shifted_copyvars(inner_to, to_start, inner_from, from_start, n, pivot)
    for i=1:n
        if i < pivot
            inner_to.x[to_start + i - 1] = inner_from.x[from_start + i - 1]
            inner_to.mult_x_L[to_start + i - 1] = inner_from.mult_x_L[from_start + i - 1]
            inner_to.mult_x_U[to_start + i - 1] = inner_from.mult_x_U[from_start + i - 1]
        else
            inner_to.x[to_start + i - 1] = inner_from.x[from_start + i]
            inner_to.mult_x_L[to_start + i - 1] = inner_from.mult_x_L[from_start + i]
            inner_to.mult_x_U[to_start + i - 1] = inner_from.mult_x_U[from_start + i]
        end
    end
end

function set_warm_from_sc(m_cur, m_sc, circuit, T; neg_g=-1)
    in_cur = internalmodel(m_cur).inner
    in_sc = internalmodel(m_sc).inner

    num_linconstr_cur = MathProgBase.numlinconstr(m_cur)
    num_linconstr_sc = MathProgBase.numlinconstr(m_sc)

    start_cur = linearindex(m_cur[:Vm][1])
    start_sc = linearindex(m_sc[:Vm][1])
    copyvars(in_cur, start_cur, in_sc, start_sc, T*length(circuit.bus))

    start_cur = linearindex(m_cur[:Va][1])
    start_sc = linearindex(m_sc[:Va][1])
    copyvars(in_cur, start_cur, in_sc, start_sc, T*length(circuit.bus))

    start_cur = linearindex(m_cur[:omega][1])
    start_sc = linearindex(m_sc[:omega][1])
    copyvars(in_cur, start_cur, in_sc, start_sc, T)

    for t=1:T
        start_cur = linearindex(m_cur[:Pg_ref][t,1])
        start_sc = linearindex(m_sc[:Pg_ref][t,1])
        shifted_copyvars(in_cur, start_cur, in_sc, start_sc, length(circuit.gen), neg_g)
    end

    for t=1:T
        start_cur = linearindex(m_cur[:Pg][t,1])
        start_sc = linearindex(m_sc[:Pg][t,1])
        shifted_copyvars(in_cur, start_cur, in_sc, start_sc, length(circuit.gen), neg_g)
    end

    for t=1:T
        start_cur = linearindex(m_cur[:Qg][t,1])
        start_sc = linearindex(m_sc[:Qg][t,1])
        shifted_copyvars(in_cur, start_cur, in_sc, start_sc, length(circuit.gen), neg_g)
    end

    start_cur = num_linconstr_cur + linearindex(m_cur[:pfreal][1])
    start_sc = num_linconstr_sc + linearindex(m_sc[:pfreal][1])
    copyconstrs(in_cur, start_cur, in_sc, start_sc, T*length(circuit.bus))

    start_cur = num_linconstr_cur + linearindex(m_cur[:pfimag][1])
    start_sc = num_linconstr_sc + linearindex(m_sc[:pfimag][1])
    copyconstrs(in_cur, start_cur, in_sc, start_sc, T*length(circuit.bus))
end

function set_warm(m_cur, m_prev, circuit, demand;
                  warm_type=:shift_copy,
                  profname="noname",
                  opt::Option = Option())
    num_buses = length(circuit.bus)
    num_gens = length(circuit.gen)
    num_linconstr = MathProgBase.numlinconstr(m_cur)

    rateA = getfield.(circuit.line, :rateA)
    num_linelimits = length(findall((rateA .!= 0) .& (rateA .< 1.0e10)))

    T = size(demand.pd, 2)
    @assert T >= 2

    inner_cur = internalmodel(m_cur).inner
    inner_prev = internalmodel(m_prev).inner

    size_pl = 0
    if opt.piecewise
        # They all have the same number of pieces so we use the first one.
        size_pl = num_gens*(circuit.gen[1].n - 1)
    end

    # -----------------------------------------------------------------
    # Copy variables and their multipliers by shifting.
    # Assume that time is in dimension 1.
    # -----------------------------------------------------------------

    start = linearindex(m_cur[:Pg][1])
    copyvars(inner_cur, start, inner_prev, start+num_gens, num_gens*(T-1))

    start = linearindex(m_cur[:Qg][1])
    copyvars(inner_cur, start, inner_prev, start+num_gens, num_gens*(T-1))

    start = linearindex(m_cur[:Vm][1])
    copyvars(inner_cur, start, inner_prev, start+num_buses, num_buses*(T-1))

    start = linearindex(m_cur[:Va][1])
    copyvars(inner_cur, start, inner_prev, start+num_buses, num_buses*(T-1))

    if opt.piecewise
        start = linearindex(m_cur[:Cg][1])
        copyvars(inner_cur, start, inner_prev, start+num_gens, num_gens*(T-1))
    end

    if opt.freq_ctrl
        start = linearindex(m_cur[:Pg_ref][1])
        copyvars(inner_cur, start, inner_prev, start+num_gens, num_gens*(T-1))

        start = linearindex(m_cur[:omega][1])
        copyvars(inner_cur, start, inner_prev, start+1, (T-1))
    end

    # -----------------------------------------------------------------
    # Copy constraint values and their multipliers by shifting.
    # Assume that time is in dimension 1.
    # -----------------------------------------------------------------

    if T > 2

        # -------------------------------------------------------------
        # Move the nonzero multipliers of the first ramping constraint
        # to the appropriate bound multipliers of variables.
        #
        # This is because the first ramping constraint was merged into
        # the bound constraints as we shift the time horizon.
        #
        #   e.g., setlower(Pg[2], max(Pg[2].Pmin, Pg[1] - ramping))
        #         setupper(Pg[2], min(Pg[2].Pmax, Pg[1] + ramping))
        #
        # The sign of a multiplier depends on the direction of its
        # constraint:
        #
        #  '>=' has a nonpositive sign.
        #  '<=' has a nonnegative sign.
        #  'a <= g(x) <= b' is equivalent to g(x) >= a and g(x) <= b.
        #   Then the sign follows the previous rule.
        #
        # Also a constraint value depends on the type.
        #
        #   nonlinear constraint: value of g == g(x) - rhs if one-sided,
        #                         value of g == g(x) if double-sided.
        #   linear constraint   : value of g == g(x)
        #
        # For multipliers on variable bounds, they all have nonnegative
        # signs.
        #
        # There is no consistency in JuMP.
        # -------------------------------------------------------------

        if opt.freq_ctrl
            Pgstart = linearindex(m_cur[:Pg_ref][1])
        else
            Pgstart = linearindex(m_cur[:Pg][1])
        end
        tol = 1e-6

        start = linearindex(m_cur[:ramping][1])
        copyconstrs(inner_cur, start, inner_prev,
                    start+num_gens, num_gens*(T-2))

        # ---------------------------------------------------------
        # Identify active ramping constraints and copy their
        # multipliers to the appropriate bound multipliers.
        # ---------------------------------------------------------

        for g=1:num_gens
            val_g = inner_prev.g[start+g-1]
            mult_g = inner_prev.mult_g[start+g-1]
            ramp_agc = circuit.gen[g].ramp_agc

            if val_g <=  -ramp_agc + tol && mult_g < -tol

                # -------------------------------------------------
                # Ramping down has occured, thus lower bound is
                # active. Copy the multiplier to lower bound mult.
                # -------------------------------------------------

                inner_cur.mult_x_L[Pgstart+g-1] = -mult_g
            elseif val_g >= ramp_agc - tol && mult_g > tol

                # -------------------------------------------------
                # Ramping up has occured, thus upper bound is
                # active. Copy the multiplier to upper bound mult.
                # -------------------------------------------------

                inner_cur.mult_x_U[Pgstart+g-1] = mult_g
            end
        end
    end

    if opt.piecewise
        start = linearindex(m_cur[:plcurve][1,1,1])
        copyconstrs(inner_cur, start, inner_prev, start+size_pl, size_pl*(T-1))
    end

    if opt.freq_ctrl
        start = linearindex(m_cur[:pg_total][1])
        copyconstrs(inner_cur, start, inner_prev, start+num_gens, num_gens*(T-1))
    end

    start = num_linconstr + linearindex(m_cur[:pfreal][1])
    copyconstrs(inner_cur, start, inner_prev, start+num_buses, num_buses*(T-1))

    start = num_linconstr + linearindex(m_cur[:pfimag][1])
    copyconstrs(inner_cur, start, inner_prev, start+num_buses, num_buses*(T-1))

    if num_linelimits > 0
        start = num_linconstr + linearindex(m_cur[:flowmaxfrom][1])
        copyconstrs(inner_cur, start, inner_prev,
                    start+num_linelimits, num_linelimits*(T-1))

        start = num_linconstr + linearindex(m_cur[:flowmaxto][1])
        copyconstrs(inner_cur, start, inner_prev,
                    start+num_linelimits, num_linelimits*(T-1))
    end

    if warm_type == :shift_copy

        # -----------------------------------------------------------------
        # Copy the (T-1)th values for the incoming period.
        # -----------------------------------------------------------------

        start = linearindex(m_cur[:Pg][T,1])
        copyvars(inner_cur, start, inner_prev, start, num_gens)

        start = linearindex(m_cur[:Qg][T,1])
        copyvars(inner_cur, start, inner_prev, start, num_gens)

        start = linearindex(m_cur[:Vm][T,1])
        copyvars(inner_cur, start, inner_prev, start, num_buses)

        start = linearindex(m_cur[:Va][T,1])
        copyvars(inner_cur, start, inner_prev, start, num_buses)

        if opt.piecewise
            start = linearindex(m_cur[:Cg][T,1])
            copyvars(inner_cur, start, inner_prev, start, num_gens)

            start = linearindex(m_cur[:plcurve][T,1,1])
            copyconstrs(inner_cur, start, inner_prev, start, size_pl)
        end

        if opt.freq_ctrl
            start = linearindex(m_cur[:Pg_ref][T,1])
            copyvars(inner_cur, start, inner_prev, start, num_gens)

            start = linearindex(m_cur[:omega][T])
            copyvars(inner_cur, start, inner_prev, start, 1)

            start = linearindex(m_cur[:pg_total][T,1])
            copyconstrs(inner_cur, start, inner_prev, start, num_gens)
        end

        start = linearindex(m_cur[:ramping][T-1,1])
        copyconstrs(inner_cur, start, inner_prev, start, num_gens)

        # -----------------------------------------------------------------
        # Compute the initial KKT error.
        # -----------------------------------------------------------------

        d_err = max(maximum(abs.(demand.pd[:,T] .- demand.pd[:,T-1])),
                    maximum(abs.(demand.qd[:,T] .- demand.qd[:,T-1])))

        ω_err = 0
        for i=1:num_gens
            if ω_err < abs(inner_prev.mult_g[start + i - 1])
                ω_err = abs(inner_prev.mult_g[start + i - 1])
            end
        end

        println("---------------------------------------------------------")
        println("Initial error: |Δd| = ", d_err, " |ω| = ", ω_err)
        println("---------------------------------------------------------")
        flush(stdout)

        start = num_linconstr + linearindex(m_cur[:pfreal][T,1])
        copyconstrs(inner_cur, start, inner_prev, start, num_buses)

        start = num_linconstr + linearindex(m_cur[:pfimag][T,1])
        copyconstrs(inner_cur, start, inner_prev, start, num_buses)

        if num_linelimits > 0
            start = num_linconstr + linearindex(m_cur[:flowmaxfrom][T,1])
            copyconstrs(inner_cur, start, inner_prev, start, num_linelimits)

            start = num_linconstr + linearindex(m_cur[:flowmaxto][T,1])
            copyconstrs(inner_cur, start, inner_prev, start, num_linelimits)
        end
    elseif warm_type == :shift_phase1

        # -----------------------------------------------------------------
        # Solve a subproblem consisting of only the incoming period.
        # -----------------------------------------------------------------

        p1_demand = Load(demand.pd[:,T], demand.qd[:,T])
        start = -1
        if opt.freq_ctrl
            start = linearindex(m_prev[:Pg_ref][T,1])
        else
            start = linearindex(m_prev[:Pg][T,1])
        end
        prev_val = inner_prev.x[start:start+num_gens-1]

        p1_opt = Option()
        p1_opt.obj_gencost = opt.obj_gencost
        p1_opt.has_ramping = opt.has_ramping
        p1_opt.freq_ctrl = opt.freq_ctrl
        p1_opt.load_shed = opt.load_shed
        p1_opt.sc_constr = opt.sc_constr
        p1_opt.piecewise = opt.piecewise
        p1_opt.weight_loadshed = opt.weight_loadshed
        p1_opt.weight_freqctrl = opt.weight_freqctrl

        p1_opt.phase1 = true
        p1_opt.prev_val = prev_val

        m_p1 = get_mpcmodel(circuit, p1_demand; opt=p1_opt)

        # -----------------------------------------------------------------
        # Warm-start the subproblem.
        # -----------------------------------------------------------------

        setsolver(m_p1, IpoptSolver(option_file_name="ipopt.op2"))
        JuMP.build(m_p1)

        in_p1 = internalmodel(m_p1)
        inner_p1 = in_p1.inner
        p1_num_linconstr = MathProgBase.numlinconstr(m_p1)

        start = linearindex(m_cur[:Pg][T,1])
        p1start = linearindex(m_p1[:Pg][1])
        copyvars(inner_p1, p1start, inner_prev, start, num_gens)

        start = linearindex(m_cur[:Qg][T,1])
        p1start = linearindex(m_p1[:Qg][1])
        copyvars(inner_p1, p1start, inner_prev, start, num_gens)

        start = linearindex(m_cur[:Vm][T,1])
        p1start = linearindex(m_p1[:Vm][1])
        copyvars(inner_p1, p1start, inner_prev, start, num_buses)

        start = linearindex(m_cur[:Va][T,1])
        p1start = linearindex(m_p1[:Va][1])
        copyvars(inner_p1, p1start, inner_prev, start, num_buses)

        if opt.piecewise
            start = linearindex(m_cur[:Cg][T,1])
            p1start = linearindex(m_p1[:Cg][1])
            copyvars(inner_p1, p1start, inner_prev, start, num_gens)

            start = linearindex(m_cur[:plcurve][T,1,1])
            p1start = linearindex(m_p1[:plcurve][1,1,1])
            copyconstrs(inner_p1, p1start, inner_prev, start, size_pl)
        end

        if opt.freq_ctrl
            start = linearindex(m_cur[:Pg_ref][T,1])
            p1start = linearindex(m_p1[:Pg_ref][1])
            copyvars(inner_p1, p1start, inner_prev, start, num_gens)

            start = linearindex(m_cur[:omega][T])
            p1start = linearindex(m_p1[:omega][1])
            copyvars(inner_p1, p1start, inner_prev, start, 1)

            start = linearindex(m_cur[:pg_total][T,1])
            p1start = linearindex(m_p1[:pg_total][1])
            copyconstrs(inner_p1, p1start, inner_prev, start, num_gens)
        end

        start = linearindex(m_cur[:ramping][T-1,1])
        p1start = linearindex(m_p1[:ramping][1])
        copyconstrs(inner_p1, p1start, inner_prev, start, num_gens)

        start = num_linconstr + linearindex(m_cur[:pfreal][T,1])
        p1start = p1_num_linconstr + linearindex(m_p1[:pfreal][1])
        copyconstrs(inner_p1, p1start, inner_prev, start, num_buses)

        start = num_linconstr + linearindex(m_cur[:pfimag][T,1])
        p1start = p1_num_linconstr + linearindex(m_p1[:pfimag][1])
        copyconstrs(inner_p1, p1start, inner_prev, start, num_buses)

        if num_linelimits > 0
            start = num_linconstr + linearindex(m_cur[:flowmaxfrom][T,1])
            p1start = p1_num_linconstr + linearindex(m_p1[:flowmaxfrom][1])
            copyconstrs(inner_p1, p1start, inner_prev, start, num_linelimits)

            start = num_linconstr + linearindex(m_cur[:flowmaxto][T,1])
            p1start = p1_num_linconstr + linearindex(m_p1[:flowmaxto][1])
            copyconstrs(inner_p1, p1start, inner_prev, start, num_linelimits)
        end

        MathProgBase.setwarmstart!(in_p1, inner_p1.x)
        MathProgBase.optimize!(in_p1)

        stat = MathProgBase.status(in_p1)

        if stat != :Infeasible && stat != :Unbounded
            m_p1.colVal = MathProgBase.getsolution(in_p1)
        else
            println("SPOPF has failed with status: ", stat)
            @assert false
        end

        save_rampinfo(m_p1, 2, num_gens, "rampinfo_p1_"*profname; circuit=circuit)

        # -----------------------------------------------------------------
        # Copy the values from the Phase1 solution for the incoming period.
        # -----------------------------------------------------------------

        start = linearindex(m_cur[:Pg][T,1])
        p1start = linearindex(m_p1[:Pg][1])
        copyvars(inner_cur, start, inner_p1, p1start, num_gens)

        start = linearindex(m_cur[:Qg][T,1])
        p1start = linearindex(m_p1[:Qg][1])
        copyvars(inner_cur, start, inner_p1, p1start, num_gens)

        start = linearindex(m_cur[:Vm][T,1])
        p1start = linearindex(m_p1[:Vm][1])
        copyvars(inner_cur, start, inner_p1, p1start, num_buses)

        start = linearindex(m_cur[:Va][T,1])
        p1start = linearindex(m_p1[:Va][1])
        copyvars(inner_cur, start, inner_p1, p1start, num_buses)

        if opt.piecewise
            start = linearindex(m_cur[:Cg][T,1])
            p1start = linearindex(m_p1[:Cg][1])
            copyvars(inner_cur, start, inner_p1, p1start, num_gens)

            start = linearindex(m_cur[:plcurve][T,1,1])
            p1start = linearindex(m_p1[:plcurve][1,1,1])
            copyconstrs(inner_cur, start, inner_p1, p1start, size_pl)
        end

        if opt.freq_ctrl
            start = linearindex(m_cur[:Pg_ref][T,1])
            p1start = linearindex(m_p1[:Pg_ref][1])
            copyvars(inner_cur, start, inner_p1, p1start, num_gens)

            start = linearindex(m_cur[:omega][T])
            p1start = linearindex(m_p1[:omega][1])
            copyvars(inner_cur, start, inner_p1, p1start, 1)

            start = linearindex(m_cur[:pg_total][T,1])
            p1start = linearindex(m_p1[:pg_total][1])
            copyconstrs(inner_cur, start, inner_p1, p1start, num_gens)
        end

        start = linearindex(m_cur[:ramping][T-1,1])
        p1start = linearindex(m_p1[:ramping][1])
        copyconstrs(inner_cur, start, inner_p1, p1start, num_gens)

        # -----------------------------------------------------------------
        # Compute the initial KKT error.
        # -----------------------------------------------------------------

        d_err = max(maximum(abs.(demand.pd[:,T] .- demand.pd[:,T-1])),
                    maximum(abs.(demand.qd[:,T] .- demand.qd[:,T-1])))

        ω_err = 0
        for i=1:num_gens
            if ω_err < abs(inner_p1.mult_g[p1start + i - 1])
                ω_err = abs(inner_p1.mult_g[p1start + i - 1])
            end
        end

        println("---------------------------------------------------------")
        println("Initial error: |Δd| = ", d_err, " |ω| = ", ω_err)
        println("---------------------------------------------------------")
        flush(stdout)

        start = num_linconstr + linearindex(m_cur[:pfreal][T,1])
        p1start = p1_num_linconstr + linearindex(m_p1[:pfreal][1])
        copyconstrs(inner_cur, start, inner_p1, p1start, num_buses)

        start = num_linconstr + linearindex(m_cur[:pfimag][T,1])
        p1start = p1_num_linconstr + linearindex(m_p1[:pfimag][1])
        copyconstrs(inner_cur, start, inner_p1, p1start, num_buses)

        if num_linelimits > 0
            start = num_linconstr + linearindex(m_cur[:flowmaxfrom][T,1])
            p1start = p1_num_linconstr + linearindex(m_p1[:flowmaxfrom][1])
            copyconstrs(inner_cur, start, inner_p1, p1start, num_linelimits)

            start = num_linconstr + linearindex(m_cur[:flowmaxto][T,1])
            p1start = p1_num_linconstr + linearindex(m_p1[:flowmaxto][1])
            copyconstrs(inner_cur, start, inner_p1, p1start, num_linelimits)
        end
    end
end

function save_header(m, circuit, T, basename="solution_header")
    num_linconstr = MathProgBase.numlinconstr(m)
    num_buses = length(circuit.bus)
    num_gens = length(circuit.gen)
    rateA = getfield.(circuit.line, :rateA)
    num_linelimits = length(findall((rateA .!= 0) .& (rateA .< 1.0e10)))

    # ----------------------------------------------------------------
    # Save the relative locations.
    # ----------------------------------------------------------------

    f = open(basename*".txt", "w")
    write(f, string(num_linconstr), "\n")
    write(f, string(T), "\n")
    write(f, string(num_gens), "\n")
    write(f, string(num_buses), "\n")
    write(f, string(num_linelimits), "\n")

    start = linearindex(m[:Pg][1])
    write(f, string(start), "\n")

    start = linearindex(m[:Qg][1])
    write(f, string(start), "\n")

    start = linearindex(m[:Vm][1])
    write(f, string(start), "\n")

    start = linearindex(m[:Va][1])
    write(f, string(start), "\n")

    if T > 1
        start = linearindex(m[:ramping][1])
        write(f, string(start), "\n")
    end

    start = num_linconstr + linearindex(m[:pfreal][1])
    write(f, string(start), "\n")

    start = num_linconstr + linearindex(m[:pfimag][1])
    write(f, string(start), "\n")

    if num_linelimits > 0
    start = num_linconstr + linearindex(m[:flowmaxfrom][1])
    write(f, string(start), "\n")

    start = num_linconstr + linearindex(m[:flowmaxto][1])
    write(f, string(start), "\n")
    end

    close(f)
end

function save_rampagc(circuit, basename="solution_rampagc")
    gen = circuit.gen

    f = open(basename*".txt", "w")
    for g = 1:length(gen)
        write(f, string(gen[g].ramp_agc), "\n")
    end
    close(f)
end

function load_solution(m, basename="solution")
    basename = replace(basename, "cold" => "warm")
    basename = replace(basename, "phase1" => "phase0")
    basename = replace(basename, "use_qp" => "no_qp")
    basename = replace(basename, r"h[\d]+_" => "")
#    basename = replace(basename, r"_pf[\d]" => "")
    basename = replace(basename, r"pt0.[\d]+" => "pt0")
    basename = replace(basename, r"_mlim_[\d]+" => "")

    println("Loading a solution from the file: ", basename*".jld2")
    sol = load(basename*".jld2")

    if !m.internalModelLoaded
        JuMP.build(m)
    end

    @assert(sol["n"] == MathProgBase.numvar(m))
    @assert(sol["m"] == MathProgBase.numconstr(m))

    inner = internalmodel(m).inner

    copyto!(inner.x, 1, sol["x"], 1, sol["n"])
    copyto!(inner.mult_x_L, 1, sol["mult_x_L"], 1, sol["n"])
    copyto!(inner.mult_x_U, 1, sol["mult_x_U"], 1, sol["n"])
    copyto!(inner.g, 1, sol["g"], 1, sol["m"])
    copyto!(inner.mult_g, 1, sol["mult_g"], 1, sol["m"])
    m.colVal = MathProgBase.getsolution(internalmodel(m))
    m.objVal = MathProgBase.getobjval(internalmodel(m))
end

function load_circuit(basename="solution_circuit")
    basename = replace(basename, "cold" => "warm")
    basename = replace(basename, "phase1" => "phase0")
    basename = replace(basename, "use_qp" => "no_qp")
    basename = replace(basename, r"h[\d]+_" => "")
#    basename = replace(basename, r"_pf[\d]" => "")
    basename = replace(basename, r"pt0.[\d]+" => "pt0")
    basename = replace(basename, r"_mlim_[\d]+" => "")

    println("Loading a circuit from the file: ", basename*".jld2")
    dict = load(basename*".jld2")

    circuit = dict["circuit"]
    pg_old_x = dict["pg_old_x"]
    neg_g = dict["neg_g"]

    return circuit, pg_old_x, neg_g
end

function save_circuit(circuit, pg_old_x, neg_g, basename="solution_circuit")
    basename = replace(basename, r"h[\d]+_" => "")
    save(basename*".jld2", "circuit", circuit, "pg_old_x", pg_old_x, "neg_g", neg_g)
end

function save_solution(m, basename="solution")
    basename = replace(basename, r"h[\d]+_" => "")
    inner = internalmodel(m).inner
    save(basename*".jld2", "n", length(inner.x), "m", length(inner.g), "x", inner.x,
         "mult_x_L", inner.mult_x_L, "mult_x_U", inner.mult_x_U,
         "g", inner.g, "mult_g", inner.mult_g)
end

function save_rampinfo(m, T, num_gens, basename="rampinfo", mode="a";
               circuit = nothing)
    if T <= 1
        return
    end

    f = open(basename*".txt", mode)
    inner = internalmodel(m).inner

    # ----------------------------------------------------------------
    # Save the ratio of binded ramping constraints.
    # ----------------------------------------------------------------

    start = linearindex(m[:ramping][1])

    for t=1:T-1
        idx = start + (t-1)*num_gens
        num_binded = 0
        num_constr_binded = 0

        for g=1:num_gens
            if abs(inner.mult_g[idx + g - 1]) >= 1e-3
                num_binded += 1
            end

            if circuit != nothing
                if abs(abs(inner.g[idx + g - 1]) - circuit.gen[g].ramp_agc) < 1e-6
                    num_constr_binded += 1
                end
            end
        end

        s = @sprintf("(%.6f  %.6f)", num_binded/num_gens, num_constr_binded/num_gens)

        if t == 1
            write(f, s)
        else
            write(f, "\t", s)
        end
    end

    write(f, "\n")
    close(f)
end

function solve_cold(m, circuit, demand; opt = Option())
    stat = :Infeasible

    if opt.powerflow_solve
        m_pf = get_mpcpfmodel(circuit, demand)
        setsolver(m_pf, IpoptSolver(option_file_name="ipopt.opt"))
        stat = solve(m_pf)
        if stat != :Optimal
            println("Power flow stat is not optimal: ", stat)
            return
        end

        setsolver(m, IpoptSolver(option_file_name="ipopt.pf"))
        JuMP.build(m)
        MathProgBase.setwarmstart!(internalmodel(m), m_pf.colVal)
        MathProgBase.optimize!(internalmodel(m))
        stat = MathProgBase.status(internalmodel(m))
        if stat != :Infeasible && stat != :Unbounded
            T = size(demand.pd, 2)
            m.colVal = MathProgBase.getsolution(internalmodel(m))
            m.objVal = get_mpcobjectivevalue(m, circuit, T)
        end
    else
        stat = solve(m)
    end

    return stat
end

function print_statistics(m, circuit, T, opt)
    baseMVA = circuit.baseMVA
    num_buses = length(circuit.bus)
    num_gens = length(circuit.gen)

    if opt.freq_ctrl
        println("\nOmega and Frequencies:")
        for t in 1:T
            @printf("%.16f    %.16f\n",
                    getvalue(m[:omega][t]), 60 + 60*getvalue(m[:omega][t]))
        end

        val = 0
        for g in 1:num_gens
            val += circuit.gen[g].alpha*getvalue(m[:omega][1])
        end

        @printf("\nSum over g of alpha*omega: %.6e\n", val)
    end

    if opt.load_shed
        @printf("\nLoad %10s %10s:\n", "Shedding", "Surplus")
        for t=1:T
            shed = 0
            surp = 0

            for b in 1:num_buses
                shed += getvalue(m[:sigma_p][t, b])*baseMVA
                surp += getvalue(m[:sigma_m][t, b])*baseMVA
            end

            @printf("     %.6e %.6e\n", shed, surp)
        end
    end

    gen_cost = mpc_get_gencost(m, circuit, T; opt = opt)
    freq_cost = mpc_get_freqcost(m, circuit, T; opt = opt)
    shed_cost = mpc_get_shedcost(m, circuit, T; opt = opt)

    @printf("Objective value:\n")
    @printf("  Generation cost = %.6e\n", gen_cost)
    @printf("  Frequency  cost = %.6e\n", freq_cost)
    @printf("  Shedding   cost = %.6e\n", shed_cost)
    @printf("  Total      cost = %.6e\n", getobjectivevalue(m))

    return
end

function usage()
    println("Usage: julia mpc.jl case scen T H LS RS warm opt [profname]")
    println(" where")
    println("          case - the name of the case file")
    println("          scen - the name of the scenario file")
    println("             T - the length of the time horizon")
    println("             H - the length of the time horizon shift")
    println("            LS - load scale")
    println("            RS - ramping scale")
    println("          warm - cold|shift_copy|shift_phase1")
    println("           opt - option file number for warm-start")
    println("      cut_line - 1 to cut a line, 0 otherwise")
    println("       cut_gen - 1 to turn off a generator, 0 otherwise")
    println("       perturb - perturbation scale")
    println("            qp - 1 if qp approx is used 0 otherwise")
    println("      load_sol - 1 if solution for time horizon 1 is available")
    println("      pf_solve - 1 if want to solve power flow for init point")
    println("      profname - profile name")
end

function main(args)

    # ---------------------------------------------------------------------
    # Parse the arguments.
    # ---------------------------------------------------------------------

    if length(args) < 14
        usage()
        return
    end

    case = args[1]
    scen = args[2]

    T = max(parse(Int,args[3]),1)
    H = max(parse(Int,args[4]),1)
    load_scale = parse(Float64,args[5])
    ramp_scale = parse(Float64,args[6])
    warm = args[7]
    optfile_num = max(parse(Int,args[8]),2)
    cut_line = (parse(Int,args[9]) == 1) ? true : false
    cut_gen = (parse(Int,args[10]) == 1) ? true : false
    perturb = parse(Float64,args[11])
    qp = parse(Int,args[12])
    load_sol = parse(Int,args[13])
    powerflow_solve = (parse(Int,args[14]) == 1) ? true : false
    profname = nothing

    if length(args) == 15
        profname = args[15]
    end

    warmstart = true
    if warm == "cold"
        warmstart = false
    elseif warm == "shift_copy"
        warm_type = :shift_copy
    elseif warm == "shift_phase1"
        warm_type = :shift_phase1
    else
        usage()
        return
    end

    if cut_line == 1 && cut_gen == 1
        println("Error: a line and a generator cannot be both turned off.")
        usage()
        return
    end

    baseMVA = 100
    piecewise_cost = false

    println("Options specified:")
    println("        case: ", case)
    println("        scen: ", scen)
    println("           T: ", T)
    println("           H: ", H)
    println("  load scale: ", load_scale)
    println("  ramp scale: ", ramp_scale)
    println("  warm-start: ", warm)
    println("         opt: ", optfile_num)
    println("    cut_line: ", cut_line)
    println("     cut_gen: ", cut_gen)
    println("     perturb: ", perturb)
    println("          qp: ", qp)
    println("    load_sol: ", load_sol)
    println("    pf_solve: ", powerflow_solve)
    println("    profname: ", (profname != nothing) ? profname : "default")
    flush(stdout)

    if profname == nothing
        profname = get_profname(case, T, load_scale, ramp_scale, warm)
    end

    # ---------------------------------------------------------------------
    # Read the circuit and load.
    # ---------------------------------------------------------------------

    circuit = getcircuit(case, baseMVA, ramp_scale)
    load = getload(scen, load_scale)

    println("Network statistics:")
    println("   # buses     : ", length(circuit.bus))
    println("   # generators: ", length(circuit.gen))
    println("   # branches  : ", length(circuit.line))

    if length(circuit.gen) > 0 && circuit.gen[1].gentype == 1
        piecewise_cost = true
    end

    println("Number of lines: ", length(circuit.line))
    flush(stdout)

    # ---------------------------------------------------------------------
    # If cut_line is set, we solve the problem and cut a line having
    # the least flow. Then re-read the circuit.
    # ---------------------------------------------------------------------

    @assert(cut_line || cut_gen)
    pg_old_x = Float64[]
    weight_freqctrl = 1e3

    if load_sol != 1
        neg_line = neg_gen = -1

        opt = Option()
        opt.obj_gencost = true
        opt.has_ramping = false
        opt.piecewise = piecewise_cost
        demand = Load(load.pd[:,1], load.qd[:,1])
        m_cur = get_mpcmodel(circuit, demand; opt=opt)
        init_x(m_cur, circuit, demand)
        setsolver(m_cur, IpoptSolver(option_file_name="ipopt.opt"))
        stat = solve(m_cur)

        if stat != :Optimal
            println("Stat is not optimal: ", stat)
            return
        end

        if cut_line
            neg_line = get_cutline(m_cur, circuit, 1)

            # Manually choose a non-zero branch to cut because there are cases
            # where the smallest one has zero value:
            #    6 (5019,9112) for 1354pegase
            #   35 ( 346,  10) for 2383wp
            #    3 (2169,3124) for 9241pegase
            if occursin("case9/", scen)
                neg_line = 2
            elseif occursin("case30/", scen)
                neg_line = 12
            elseif occursin("case300/", scen)
                neg_line = 5
            elseif occursin("1354pegase", case)
                neg_line = 6
            elseif occursin("2383wp", case)
                neg_line = 35
            elseif occursin("9241pegase", case)
                neg_line = 3
            end

            @printf("Cut line (%5d,%5d) . . .\n",
                    circuit.line[neg_line].from, circuit.line[neg_line].to)
        end

        if cut_gen
            neg_gen = get_cutgen(m_cur, circuit, 1)
            @printf("Turn off generator %d . . .\n", neg_gen)
        end
        flush(stdout)

        opt = Option()
        opt.obj_gencost = true
        opt.has_ramping = true
        opt.freq_ctrl = true
        opt.sc_constr = true
        opt.neg_g = [neg_gen]
        opt.weight_freqctrl = weight_freqctrl

        demand = Load(load.pd[:,1:T], load.qd[:,1:T])
        m_sc = get_mpcmodel(circuit, demand; opt = opt)
        init_x(m_sc, circuit, demand; opt = opt)
        setsolver(m_sc, IpoptSolver(option_file_name="ipopt.opt"))
        stat = solve(m_sc)
        if stat != :Optimal
            println("Stat is not optimal: ", stat)
            return
        end

        println("\nOmega and Frequencies:")
        for t in 1:T
            @printf("%.6e    %.16f\n",
                    getvalue(m_sc[:omega][t]), 60 + 60*getvalue(m_sc[:omega][t]))
        end

        pos_gen = findall(x -> x != neg_gen, collect(1:length(circuit.gen)))
        pg_old_x = zeros(length(pos_gen))
        for (j,g) in enumerate(pos_gen)
           pg_old_x[j] = getvalue(m_sc[:Pg_f][g])
        end

        circuit = getcircuit(case, baseMVA, ramp_scale;
                             neg_line=neg_line, neg_gen=neg_gen)
        save_circuit(circuit, pg_old_x, [neg_gen], "solution_circuit_"*profname)
    else
        circuit, pg_old_x, neg_g = load_circuit("solution_circuit_"*profname)
    end

    println("Number of lines: ", length(circuit.line))
    flush(stdout)

    num_gens = length(circuit.gen)
    num_buses = length(circuit.bus)
    gen = circuit.gen

    # ---------------------------------------------------------------------
    # Solve the first time horizon [1:T] using Ipopt.
    # ---------------------------------------------------------------------

    single_objval = 0
    total_objval = 0
    total_single_objval = 0

    opt = Option()
    opt.obj_gencost = false
    opt.freq_ctrl = true
    opt.sc_constr = false
    opt.piecewise = piecewise_cost
    opt.powerflow_solve = powerflow_solve
    opt.neg_g = Int[]
    opt.weight_freqctrl = weight_freqctrl

    demand = Load(load.pd[:,1:T], load.qd[:,1:T])
    m_cur = get_mpcmodel(circuit, demand; opt=opt)
    setsolver(m_cur, IpoptSolver(option_file_name="ipopt.opt"))

    if load_sol == 1
        stat = :Optimal
        load_solution(m_cur, "solution_"*profname)
        save_rampinfo(m_cur, T, num_gens, "rampinfo_"*profname, "w"; circuit=circuit)
    else
        if opt.freq_ctrl
            for g in 1:num_gens
                setlowerbound(m_cur[:Pg_ref][1,g], pg_old_x[g])
                setupperbound(m_cur[:Pg_ref][1,g], pg_old_x[g])
                setvalue(m_cur[:Pg_ref][1,g], pg_old_x[g])
            end
        else
            for g in 1:num_gens
                setlowerbound(m_cur[:Pg][1,g], pg_old_x[g])
                setupperbound(m_cur[:Pg][1,g], pg_old_x[g])
                setvalue(m_cur[:Pg][1,g], pg_old_x[g])
            end
        end

        JuMP.build(m_cur)
        set_warm_from_sc(m_cur, m_sc, circuit, T; neg_g=neg_gen)
        in_cur = internalmodel(m_cur)
        MathProgBase.setwarmstart!(in_cur, in_cur.inner.x)
        MathProgBase.optimize!(in_cur)

        stat = MathProgBase.status(in_cur)
        if stat != :Optimal
            println("Stat is not optimal: ", stat)
            return
        end

        m_cur.objVal = MathProgBase.getobjval(in_cur)
        m_cur.colVal = MathProgBase.getsolution(in_cur)

        save_solution(m_cur, "solution_"*profname)
        save_rampinfo(m_cur, T, num_gens, "rampinfo_"*profname, "w"; circuit=circuit)
    end

    single_objval = get_mpcobjectivevalue(m_cur, circuit, 1)
    total_single_objval += single_objval
    total_objval += single_objval

    print_statistics(m_cur, circuit, T, opt)

    # ---------------------------------------------------------------------
    # Warm-start option for Ipopt.
    # ---------------------------------------------------------------------

    m_prev = nothing
    for h in 2:H
        println("\nh = h", h, "\n")
        flush(stdout)

        m_prev = m_cur
        inner_prev = internalmodel(m_prev).inner

        # -----------------------------------------------------------------
        # Build a model for the current horizon: [h,h+T-1]
        # -----------------------------------------------------------------

        demand = Load(load.pd[:,h:h+T-1], load.qd[:,h:h+T-1])
        if perturb != 0
            Random.seed!(0)
            factor = rand((-1,1))*(rand(Float64,1)*perturb)
            demand.pd .= demand.pd .+ factor.*demand.pd
            demand.qd .= demand.qd .+ factor.*demand.qd
        end

        m_cur = get_mpcmodel(circuit, demand; opt=opt)

        # -----------------------------------------------------------------
        # Change the lower/upper bounds on the first variable based on
        # the solution of the previous model.
        # -----------------------------------------------------------------

        Pg_prev = Pg_cur = nothing

        if opt.freq_ctrl
            Pg_prev = m_prev[:Pg_ref]
            Pg_cur = m_cur[:Pg_ref]
        else
            Pg_prev = m_prev[:Pg]
            Pg_cur = m_cur[:Pg]
        end

        start = linearindex(Pg_prev[1])
        for g=1:num_gens
            setlowerbound(Pg_cur[1,g],
                          max(gen[g].Pmin, inner_prev.x[start+g-1] - gen[g].ramp_agc))
            setupperbound(Pg_cur[1,g],
                          min(gen[g].Pmax, inner_prev.x[start+g-1] + gen[g].ramp_agc))
        end

        if T >= 2 && warmstart
            if h == 2
                rm("rampinfo_p1_"*profname*".txt", force=true)
            end

            # -------------------------------------------------------------
            # For a warm-start,
            #
            #  - set options for Ipopt.
            #  - build a model.
            #  - copy primal/dual values and constraints values.
            # -------------------------------------------------------------

            optname = "ipopt.op"*string(optfile_num)
            setsolver(m_cur, IpoptSolver(option_file_name=optname))
            JuMP.build(m_cur)
            set_warm(m_cur, m_prev, circuit, demand;
                     warm_type=warm_type, profname=profname, opt=opt)

            if qp == 1
                d, m_qp = mpc_qp(m_cur; optfile=optname)
                copyto!(m_cur.colVal, 1, m_qp.colVal, 1, MathProgBase.numvar(m_qp))
                m_cur.objVal = get_mpcobjectivevalue(m_cur, circuit, T; opt=opt)
                constr_viol = get_mpcconstrviolation(m_cur, d)

                @printf("QP approximation objective..................: %18.16e\n", m_cur.objVal)
                @printf("QP approximation constraint violation.......: %18.16e\n", constr_viol)

                copyvars(internalmodel(m_cur).inner, 1,
                         internalmodel(m_qp).inner, 1,
                         MathProgBase.numvar(m_qp))
                copyconstrs(internalmodel(m_cur).inner, 1,
                            internalmodel(m_qp).inner, 1,
                            MathProgBase.numconstr(m_qp))
            else
                # -------------------------------------------------------------
                # Call optimize!() for a warm-start instead of solve() because
                # solve() will build a new internal model that ignores the
                # multipliers.
                #
                # Note that MathProgBase.optimize!() doesn't change the original
                # model's solution status, such as objVal, colVal, and redCosts.
                # If we need to access them using the original model, we should
                # manually assign them.
                # -------------------------------------------------------------

                in_cur = internalmodel(m_cur)   # IpoptMathProgModel is returned.
                MathProgBase.setwarmstart!(in_cur, in_cur.inner.x)
                MathProgBase.optimize!(in_cur)
                stat = MathProgBase.status(in_cur)

                if stat != :Infeasible && stat != :Unbounded
                    m_cur.objVal = MathProgBase.getobjval(in_cur)
                    m_cur.colVal = MathProgBase.getsolution(in_cur)
                end
            end
        else

            # -------------------------------------------------------------
            # Cold-start
            # -------------------------------------------------------------

            stat = solve_cold(m_cur, circuit, demand; powerflow_solve = powerflow_solve)
        end

        if stat != :Optimal && stat != :UserLimit
            println("Stat is not optimal: ", stat)
            return
        end

        save_rampinfo(m_cur, T, num_gens, "rampinfo_"*profname; circuit=circuit)

        single_objval = get_mpcobjectivevalue(m_cur, circuit, 1)
        total_single_objval += single_objval
        total_objval += single_objval

        print_statistics(m_cur, circuit, T, opt)
    end

    # Add the remaining cost.
    total_objval += getobjectivevalue(m_cur) - single_objval

    @printf("Total objective value..............: %18.6e\n", total_objval)
    @printf("Total single objective value.......: %18.6e\n", total_single_objval)
    return
end

main(ARGS)

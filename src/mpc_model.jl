function get_mpcmodel(circuit, demand;
                      opt::Option=Option())

    m = Model()

    # Shortcuts
    baseMVA = circuit.baseMVA
    busref = circuit.busref
    bus = circuit.bus
    line = circuit.line
    gen = circuit.gen
    yline = circuit.yline
    ybus = circuit.ybus
    busdict = circuit.busdict
    frombus = circuit.frombus
    tobus = circuit.tobus
    bus2gen = circuit.bus2gen

    Pd = demand.pd
    Qd = demand.qd
    T = size(Pd,2)

    num_buses = length(bus)
    num_gens = length(gen)
    num_lines = length(line)

    @variable(m, gen[g].Pmin <= Pg[t=1:T,g=1:num_gens] <= gen[g].Pmax)
    @variable(m, gen[g].Qmin <= Qg[t=1:T,g=1:num_gens] <= gen[g].Qmax)
    @variable(m, bus[b].Vmin <= Vm[t=1:T,b=1:num_buses] <= bus[b].Vmax)
    @variable(m, Va[t=1:T,b=1:num_buses])

    active_gen = collect(1:num_gens)
    if !isempty(opt.neg_g)
        active_gen = findall(x -> !(x in opt.neg_g), collect(1:num_gens))
    end
    num_active_gen = length(active_gen)

    if opt.freq_ctrl
        # omega for primary control (frequency deviation in p.u.)
        # Pg_ref for secondary control (speed-changer setting)
        # sigma_p and sigma_m (load shedding and surplus)
        @variable(m, gen[g].Pmin <= Pg_ref[t=1:T, g=1:num_gens] <= gen[g].Pmax)
        @variable(m, -1 <= omega[t=1:T] <= 1, start=0)
    end

    if opt.load_shed
        @variable(m, sigma_p[t=1:T, b=1:num_buses] >= 0, start=0)
        @variable(m, sigma_m[t=1:T, b=1:num_buses] >= 0, start=0)
    end

    for t in 1:T
        set_lower_bound(Va[t,busref], bus[busref].Va)
        set_upper_bound(Va[t,busref], bus[busref].Va)

        # setlowerbound(Va[t,busref], bus[busref].Va)
        # setupperbound(Va[t,busref], bus[busref].Va)
    end

    if !isempty(opt.neg_g)
        for t=1:T, g in opt.neg_g
            set_lower_bound(Pg[t,g], 0)
            set_upper_bound(Pg[t,g], 0)
            set_value(Pg[t,g], 0)
            set_lower_bound(Qg[t,g], 0)
            set_upper_bound(Qg[t,g], 0)
            set_value(Qg[t,g], 0)

            if opt.freq_ctrl
                set_lower_bound(Pg_ref[t,g], 0)
                set_upper_bound(Pg_ref[t,g], 0)
                set_value(Pg_ref[t,g], 0)
            end
        end
    end

    if opt.sc_constr
        @variable(m, gen[g].Pmin <= Pg_f[g=1:num_gens] <= gen[g].Pmax)
        @variable(m, gen[g].Qmin <= Qg_f[g=1:num_gens] <= gen[g].Qmax)
        @variable(m, bus[b].Vmin <= Vm_f[b=1:num_buses] <= bus[b].Vmax)
        @variable(m, Va_f[b=1:num_buses])

        set_lower_bound(Va_f[busref], bus[busref].Va)
        set_upper_bound(Va_f[busref], bus[busref].Va)

        if opt.freq_ctrl
            @constraint(m, pg_fix_first[g=1:num_active_gen], Pg_ref[1,active_gen[g]] == Pg_f[active_gen[g]])
        else
            @constraint(m, pg_fix_first[g=1:num_active_gen], Pg[1,active_gen[g]] == Pg_f[active_gen[g]])
        end
    end

    if opt.obj_gencost
        if opt.sc_constr
            if opt.piecewise
                @variable(m, Cg[g=1:num_gens])
                @constraint(m, plcurve[g=1:num_gens,p=1:gen[g].n-1],
                            Cg[g] - (((gen[g].coeff[2*p+2] - gen[g].coeff[2*p])/(gen[g].coeff[2*p+1] - gen[g].coeff[2*p-1]))*(baseMVA*Pg_f[g] - gen[g].coeff[2*p-1]) + gen[g].coeff[2*p]) >= 0
                           )
                @NLexpression(m, obj_gencost, sum(Cg[g] for g=1:num_gens))
            else
                @NLexpression(m, obj_gencost,
                              sum(gen[g].coeff[gen[g].n-2]*(baseMVA*Pg_f[g])^2
                                + gen[g].coeff[gen[g].n-1]*(baseMVA*Pg_f[g])
                                + gen[g].coeff[gen[g].n] for g=1:num_gens))
            end
        else
            if opt.piecewise
                @variable(m, Cg[t=1:T,g=1:num_gens])
                @constraint(m, plcurve[t=1:T,g=1:num_gens,p=1:gen[g].n-1],
                    Cg[t,g] - (((gen[g].coeff[2*p+2] - gen[g].coeff[2*p])/(gen[g].coeff[2*p+1] - gen[g].coeff[2*p-1]))*(baseMVA*Pg[t,g] - gen[g].coeff[2*p-1]) + gen[g].coeff[2*p]) >= 0
                           )
                @NLexpression(m, obj_gencost, sum(Cg[t,g] for t=1:T,g=1:num_gens))
            else
                @NLexpression(m, obj_gencost,
                              sum(gen[g].coeff[gen[g].n-2]*(baseMVA*Pg[t,g])^2
                                + gen[g].coeff[gen[g].n-1]*(baseMVA*Pg[t,g])
                                + gen[g].coeff[gen[g].n] for t=1:T,g=1:num_gens))
            end
        end
    else
        @NLexpression(m, obj_gencost, 0)
    end

    if opt.freq_ctrl
        @NLexpression(m, obj_freq_ctrl, 0.5*sum(omega[t]^2 for t=1:T))
    else
        @NLexpression(m, obj_freq_ctrl, 0)
    end

    if opt.load_shed
        @NLexpression(m, obj_load_shed,
                      sum(baseMVA*sigma_p[t, b] +
                          baseMVA*sigma_m[t, b] for t=1:T,b=1:num_buses))
    else
        @NLexpression(m, obj_load_shed, 0)
    end

    @NLobjective(m, Min,
                 obj_gencost + opt.weight_loadshed*obj_load_shed +
                 opt.weight_freqctrl*obj_freq_ctrl)

    if opt.freq_ctrl
        @constraint(m, pg_total[t=1:T,g=1:num_active_gen],
                       Pg[t,active_gen[g]] == Pg_ref[t,active_gen[g]] + gen[active_gen[g]].alpha*omega[t])
    end

    # Ramping up/down constraints
    if opt.has_ramping
        Pg_ramp = m[:Pg]
        if opt.freq_ctrl
            Pg_ramp = m[:Pg_ref]
        end

        if opt.phase1 == false
            @constraint(m, ramping[t=1:T-1,g=1:num_active_gen],
                        -gen[active_gen[g]].ramp_agc <= Pg_ramp[t+1,active_gen[g]] - Pg_ramp[t,active_gen[g]] <= gen[active_gen[g]].ramp_agc)
        else
            @constraint(m, ramping[t=1:T,g=1:num_active_gen],
                        -gen[active_gen[g]].ramp_agc <= Pg_ramp[t,active_gen[g]] - opt.prev_val[active_gen[g]] <= gen[active_gen[g]].ramp_agc)
        end
    end

    if opt.load_shed
        @NLexpression(m, shed[t=1:T, b=1:num_buses], sigma_p[t, b])
        @NLexpression(m, surp[t=1:T, b=1:num_buses], sigma_m[t, b])
    else
        @NLexpression(m, shed[t=1:T, b=1:num_buses], 0)
        @NLexpression(m, surp[t=1:T, b=1:num_buses], 0)
    end

    if opt.freq_ctrl
        #@NLexpression(m, lc[t=1:T, b=1:num_buses], baseMVA*bus[b].D*omega[t])
        @NLexpression(m, lc[t=1:T, b=1:num_buses], 0)
    else
        @NLexpression(m, lc[t=1:T, b=1:num_buses], 0)
    end

    # Power flow constraints: real part

    @NLconstraint(m, pfreal[t=1:T,b=1:num_buses],
                  (sum(yline[l].YffR for l in frombus[b])
                   + sum(yline[l].YttR for l in tobus[b])
                   + ybus[b].YshR)*Vm[t,b]^2
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].to]]*
                        (yline[l].YftR*cos(Va[t,b]-Va[t,busdict[line[l].to]])
                         + yline[l].YftI*sin(Va[t,b]-Va[t,busdict[line[l].to]]))
                        for l in frombus[b])
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].from]]*
                        (yline[l].YtfR*cos(Va[t,b]-Va[t,busdict[line[l].from]])
                         + yline[l].YtfI*sin(Va[t,b]-Va[t,busdict[line[l].from]]))
                        for l in tobus[b])
                  - (sum(baseMVA*Pg[t,g] for g in bus2gen[b]) - (Pd[b,t] + lc[t,b])) / baseMVA
                  == shed[t, b] - surp[t, b])

    # Power flow constraints: imaginary part
    @NLconstraint(m, pfimag[t=1:T,b=1:num_buses],
                  (sum(-yline[l].YffI for l in frombus[b])
                   + sum(-yline[l].YttI for l in tobus[b])
                   - ybus[b].YshI)*Vm[t,b]^2
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].to]]*
                        (-yline[l].YftI*cos(Va[t,b]-Va[t,busdict[line[l].to]])
                         + yline[l].YftR*sin(Va[t,b]-Va[t,busdict[line[l].to]]))
                        for l in frombus[b])
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].from]]*
                        (-yline[l].YtfI*cos(Va[t,b]-Va[t,busdict[line[l].from]])
                         + yline[l].YtfR*sin(Va[t,b]-Va[t,busdict[line[l].from]]))
                        for l in tobus[b])
                  - (sum(baseMVA*Qg[t,g] for g in bus2gen[b]) - Qd[b,t]) / baseMVA
                  == 0)

    # Line limits
    rateA = getfield.(line, :rateA)
    limind = findall((rateA .!= 0) .& (rateA .< 1.0e10))
    num_linelimits = length(limind)

    Yff_abs2 = zeros(num_linelimits)
    Yft_abs2 = zeros(num_linelimits)
    Yre = zeros(num_linelimits)
    Yim = zeros(num_linelimits)
    flowmax = zeros(num_linelimits)

    for i in 1:num_linelimits
        # Apparent power limits (from bus)
        l = limind[i]
        flowmax[i] = (line[l].rateA / baseMVA)^2
        Yff_abs2[i] = yline[l].YffR^2 + yline[l].YffI^2
        Yft_abs2[i] = yline[l].YftR^2 + yline[l].YftI^2
        Yre[i] = yline[l].YffR*yline[l].YftR + yline[l].YffI*yline[l].YftI
        Yim[i] = -yline[l].YffR*yline[l].YftI + yline[l].YffI*yline[l].YftR
    end

    @NLconstraint(m, flowmaxfrom[t=1:T,i=1:num_linelimits],
                  Vm[t,busdict[line[limind[i]].from]]^2 *
                  (Yff_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
                   + Yft_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
                   + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
                   (Yre[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
                    - Yim[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
                   ) - flowmax[i] <= 0)

    Ytf_abs2 = zeros(num_linelimits)
    Ytt_abs2 = zeros(num_linelimits)

    for i in 1:num_linelimits
        # Apparent power limits (to bus)
        l = limind[i]
        Ytf_abs2[i] = yline[l].YtfR^2 + yline[l].YtfI^2
        Ytt_abs2[i] = yline[l].YttR^2 + yline[l].YttI^2
        Yre[i] = yline[l].YtfR*yline[l].YttR + yline[l].YtfI*yline[l].YttI
        Yim[i] = -yline[l].YtfR*yline[l].YttI + yline[l].YtfI*yline[l].YttR
    end

    @NLconstraint(m, flowmaxto[t=1:T,i=1:num_linelimits],
                  Vm[t,busdict[line[limind[i]].to]]^2 *
                  (Ytf_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
                   + Ytt_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
                   + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
                   (Yre[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
                    -Yim[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
                   ) - flowmax[i] <=0)

    if opt.sc_constr
        # Power flow constraints: real part

        @NLconstraint(m, pfreal_f[b=1:num_buses],
                      (sum(yline[l].YffR for l in frombus[b])
                       + sum(yline[l].YttR for l in tobus[b])
                       + ybus[b].YshR)*Vm_f[b]^2
                      + sum(Vm_f[b]*Vm_f[busdict[line[l].to]]*
                            (yline[l].YftR*cos(Va_f[b]-Va_f[busdict[line[l].to]])
                             + yline[l].YftI*sin(Va_f[b]-Va_f[busdict[line[l].to]]))
                            for l in frombus[b])
                      + sum(Vm_f[b]*Vm_f[busdict[line[l].from]]*
                            (yline[l].YtfR*cos(Va_f[b]-Va_f[busdict[line[l].from]])
                             + yline[l].YtfI*sin(Va_f[b]-Va_f[busdict[line[l].from]]))
                            for l in tobus[b])
                      - (sum(baseMVA*Pg_f[g] for g in bus2gen[b]) - (Pd[b,1])) / baseMVA
                      == 0)

        # Power flow constraints: imaginary part
        @NLconstraint(m, pfimag_f[b=1:num_buses],
                      (sum(-yline[l].YffI for l in frombus[b])
                       + sum(-yline[l].YttI for l in tobus[b])
                       - ybus[b].YshI)*Vm_f[b]^2
                      + sum(Vm_f[b]*Vm_f[busdict[line[l].to]]*
                            (-yline[l].YftI*cos(Va_f[b]-Va_f[busdict[line[l].to]])
                             + yline[l].YftR*sin(Va_f[b]-Va_f[busdict[line[l].to]]))
                            for l in frombus[b])
                      + sum(Vm_f[b]*Vm_f[busdict[line[l].from]]*
                            (-yline[l].YtfI*cos(Va_f[b]-Va_f[busdict[line[l].from]])
                             + yline[l].YtfR*sin(Va_f[b]-Va_f[busdict[line[l].from]]))
                            for l in tobus[b])
                      - (sum(baseMVA*Qg_f[g] for g in bus2gen[b]) - Qd[b,1]) / baseMVA
                      == 0)

        # Line limits
        rateA = getfield.(line, :rateA)
        limind = findall((rateA .!= 0) .& (rateA .< 1.0e10))
        num_linelimits = length(limind)

        Yff_abs2 = zeros(num_linelimits)
        Yft_abs2 = zeros(num_linelimits)
        Yre = zeros(num_linelimits)
        Yim = zeros(num_linelimits)
        flowmax = zeros(num_linelimits)

        for i in 1:num_linelimits
            # Apparent power limits (from bus)
            l = limind[i]
            flowmax[i] = (line[l].rateA / baseMVA)^2
            Yff_abs2[i] = yline[l].YffR^2 + yline[l].YffI^2
            Yft_abs2[i] = yline[l].YftR^2 + yline[l].YftI^2
            Yre[i] = yline[l].YffR*yline[l].YftR + yline[l].YffI*yline[l].YftI
            Yim[i] = -yline[l].YffR*yline[l].YftI + yline[l].YffI*yline[l].YftR
        end

        @NLconstraint(m, flowmaxfrom_f[i=1:num_linelimits],
                      Vm_f[busdict[line[limind[i]].from]]^2 *
                      (Yff_abs2[i]*Vm_f[busdict[line[limind[i]].from]]^2
                       + Yft_abs2[i]*Vm_f[busdict[line[limind[i]].to]]^2
                       + 2*Vm_f[busdict[line[limind[i]].from]]*Vm_f[busdict[line[limind[i]].to]]*
                       (Yre[i]*cos(Va_f[busdict[line[limind[i]].from]] - Va_f[busdict[line[limind[i]].to]])
                        - Yim[i]*sin(Va_f[busdict[line[limind[i]].from]] - Va_f[busdict[line[limind[i]].to]]))
                       ) - flowmax[i] <= 0)

        Ytf_abs2 = zeros(num_linelimits)
        Ytt_abs2 = zeros(num_linelimits)

        for i in 1:num_linelimits
            # Apparent power limits (to bus)
            l = limind[i]
            Ytf_abs2[i] = yline[l].YtfR^2 + yline[l].YtfI^2
            Ytt_abs2[i] = yline[l].YttR^2 + yline[l].YttI^2
            Yre[i] = yline[l].YtfR*yline[l].YttR + yline[l].YtfI*yline[l].YttI
            Yim[i] = -yline[l].YtfR*yline[l].YttI + yline[l].YtfI*yline[l].YttR
        end

        @NLconstraint(m, flowmaxto_f[i=1:num_linelimits],
                      Vm_f[busdict[line[limind[i]].to]]^2 *
                      (Ytf_abs2[i]*Vm_f[busdict[line[limind[i]].from]]^2
                       + Ytt_abs2[i]*Vm_f[busdict[line[limind[i]].to]]^2
                       + 2*Vm_f[busdict[line[limind[i]].from]]*Vm_f[busdict[line[limind[i]].to]]*
                       (Yre[i]*cos(Va_f[busdict[line[limind[i]].from]] - Va_f[busdict[line[limind[i]].to]])
                        -Yim[i]*sin(Va_f[busdict[line[limind[i]].from]] - Va_f[busdict[line[limind[i]].to]]))
                       ) - flowmax[i] <=0)
    end

    return m
end

function get_mpcpfmodel(circuit, demand)
    m = Model()

    # Shortcuts
    baseMVA = circuit.baseMVA
    busref = circuit.busref
    bus = circuit.bus
    line = circuit.line
    gen = circuit.gen
    yline = circuit.yline
    ybus = circuit.ybus
    busdict = circuit.busdict
    frombus = circuit.frombus
    tobus = circuit.tobus
    bus2gen = circuit.bus2gen

    Pd = demand.pd
    Qd = demand.qd
    T = size(Pd,2)

    num_buses = length(bus)
    num_gens = length(gen)
    num_lines = length(line)

    @variable(m, gen[g].Pmin <= Pg[t=1:T,g=1:num_gens] <= gen[g].Pmax)
    @variable(m, gen[g].Qmin <= Qg[t=1:T,g=1:num_gens] <= gen[g].Qmax)
    @variable(m, bus[b].Vmin <= Vm[t=1:T,b=1:num_buses] <= bus[b].Vmax)
    @variable(m, Va[t=1:T,b=1:num_buses])

    for t in 1:T
        setlowerbound(Va[t,busref], bus[busref].Va)
        setupperbound(Va[t,busref], bus[busref].Va)
    end

    @objective(m, Min, 0)

    # Power flow constraints: real part
    @NLconstraint(m, pfreal[t=1:T,b=1:num_buses],
                  (sum(yline[l].YffR for l in frombus[b])
                   + sum(yline[l].YttR for l in tobus[b])
                   + ybus[b].YshR)*Vm[t,b]^2
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].to]]*
                        (yline[l].YftR*cos(Va[t,b]-Va[t,busdict[line[l].to]])
                         + yline[l].YftI*sin(Va[t,b]-Va[t,busdict[line[l].to]]))
                        for l in frombus[b])
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].from]]*
                        (yline[l].YtfR*cos(Va[t,b]-Va[t,busdict[line[l].from]])
                         + yline[l].YtfI*sin(Va[t,b]-Va[t,busdict[line[l].from]]))
                        for l in tobus[b])
                  - (sum(baseMVA*Pg[t,g] for g in bus2gen[b]) - Pd[b,t]) / baseMVA
                  == 0)

    # Power flow constraints: imaginary part
    @NLconstraint(m, pfimag[t=1:T,b=1:num_buses],
                  (sum(-yline[l].YffI for l in frombus[b])
                   + sum(-yline[l].YttI for l in tobus[b])
                   - ybus[b].YshI)*Vm[t,b]^2
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].to]]*
                        (-yline[l].YftI*cos(Va[t,b]-Va[t,busdict[line[l].to]])
                         + yline[l].YftR*sin(Va[t,b]-Va[t,busdict[line[l].to]]))
                        for l in frombus[b])
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].from]]*
                        (-yline[l].YtfI*cos(Va[t,b]-Va[t,busdict[line[l].from]])
                         + yline[l].YtfR*sin(Va[t,b]-Va[t,busdict[line[l].from]]))
                        for l in tobus[b])
                  - (sum(baseMVA*Qg[t,g] for g in bus2gen[b]) - Qd[b,t]) / baseMVA
                  == 0)

    return m
end

function get_mpcobjectivevalue(m, circuit, T; opt = Option())
    baseMVA = circuit.baseMVA
    gen = circuit.gen
    num_gens = length(gen)
    objval = 0

    if opt.obj_gencost
        if gen[1].gentype == 1
            Cg = getvalue(m[:Cg])
            objval += sum(Cg[t,g] for t=1:T,g=1:num_gens)
        else
            Pg = getvalue(m[:Pg])
            objval += sum(gen[g].coeff[gen[g].n-2]*(baseMVA*Pg[t,g])^2 +
                        gen[g].coeff[gen[g].n-1]*(baseMVA*Pg[t,g]) +
                        gen[g].coeff[gen[g].n] for t=1:T,g=1:num_gens)
        end
    end

    if opt.freq_ctrl
        omega = getvalue(m[:omega])
        objval += opt.weight_freqctrl*(0.5*sum(omega[t]^2 for t=1:T))
    end

    return objval
end

function mpc_get_gencost(m, circuit, T; opt = Option())
    baseMVA = circuit.baseMVA
    gen = circuit.gen
    num_gens = length(gen)

    gen_cost = 0
    if opt.sc_constr
        if opt.piecewise
            gen_cost += sum(getvalue(m[:Cg][g]) for g=1:num_gens)
        else
            gen_cost += sum(gen[g].coeff[gen[g].n-2]*(baseMVA*getvalue(m[:Pg_f][g]))^2
                            + gen[g].coeff[gen[g].n-1]*(baseMVA*getvalue(m[:Pg_f][g]))
                            + gen[g].coeff[gen[g].n] for g=1:num_gens)
        end
    else
        if opt.piecewise
            gen_cost += sum(getvalue(m[:Cg][t,g]) for t=1:T,g=1:num_gens)
        else
            gen_cost += sum(gen[g].coeff[gen[g].n-2]*(baseMVA*getvalue(m[:Pg][t,g]))^2
                            + gen[g].coeff[gen[g].n-1]*(baseMVA*getvalue(m[:Pg][t,g]))
                            + gen[g].coeff[gen[g].n] for t=1:T,g=1:num_gens)
        end
    end

    return gen_cost
end

function mpc_get_freqcost(m, circuit, T; opt = Option())
    if !opt.freq_ctrl
        return 0
    end

    freq_cost = 0.5*sum(getvalue(m[:omega][t])^2 for t=1:T)
    return freq_cost
end

function mpc_get_shedcost(m, circuit, T; opt = Option())
    if !opt.load_shed
        return 0
    end

    baseMVA = circuit.baseMVA
    num_buses = length(circuit.bus)

    shed_cost = sum(baseMVA*getvalue(m[:sigma_p][t,b]) +
                    baseMVA*getvalue(m[:sigma_m][t,b])
                    for t=1:T,b=1:num_buses)
    return shed_cost
end

function get_mpcconstrviolation(m, d)
    # Take the maximum norm.
    viol = 0

    # Check constraints.
    if MathProgBase.numconstr(m) > 0
        g = zeros(MathProgBase.numconstr(m))
        MathProgBase.eval_g(d, g, m.colVal)
        g_lb, g_ub = JuMP.constraintbounds(m)

        for i=1:length(g_lb)
            if g_lb[i] != -Inf
                err = max(0, -(g[i] - g_lb[i]))
                if viol < err
                    viol = err
                end
            end

            if g_ub[i] != Inf
                err = max(0, g[i] - g_ub[i])
                if viol < err
                    viol = err
                end
            end
        end
    end

    # Check bound constraints.
    for i=1:MathProgBase.numvar(m)
        if m.colLower[i] != -Inf
            err = max(0, -(m.colVal[i] - m.colLower[i]))
            if viol < err
                viol = err
            end
        end

        if m.colUpper[i] != Inf
            err = max(0, m.colVal[i] - m.colUpper[i])
            if viol < err
                viol = err
            end
        end
    end

    return viol
end

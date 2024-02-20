using Ipopt, NLPModelsIpopt
using JuMP
using ExaModels
using Random
using DelimitedFiles, Printf
include("mpc_data.jl")
include("mpc_model.jl")
include("mpc_2.jl")

case  = "data/case1354pegase"
scen = "data/case1354pegase/halfhour_30"
T = 10
H = 1
load_scale = 1.0
ramp_scale = 0.5
warm = "cold"
optfile_num = 2

cut_line = true
cut_gen = false

perturb = 1
qp = 0
load_sol = 0
powerflow_solve = false
profname = "case9_result"

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
        em = ExaModel(m_cur)
        result = ipopt(em)
        print(result)
        # setsolver(m_cur, IpoptSolver(option_file_name="ipopt.opt"))

        # stat = solve(m_cur)

        # if stat != :Optimal
        #     println("Stat is not optimal: ", stat)
        #     return
        # end

        # if cut_line
        #     neg_line = get_cutline(m_cur, circuit, 1)

        #     # Manually choose a non-zero branch to cut because there are cases
        #     # where the smallest one has zero value:
        #     #    6 (5019,9112) for 1354pegase
        #     #   35 ( 346,  10) for 2383wp
        #     #    3 (2169,3124) for 9241pegase
        #     if occursin("case9/", scen)
        #         neg_line = 2
        #     elseif occursin("case30/", scen)
        #         neg_line = 12
        #     elseif occursin("case300/", scen)
        #         neg_line = 5
        #     elseif occursin("1354pegase", case)
        #         neg_line = 6
        #     elseif occursin("2383wp", case)
        #         neg_line = 35
        #     elseif occursin("9241pegase", case)
        #         neg_line = 3
        #     end

        #     @printf("Cut line (%5d,%5d) . . .\n",
        #             circuit.line[neg_line].from, circuit.line[neg_line].to)
        # end

        # if cut_gen
        #     neg_gen = get_cutgen(m_cur, circuit, 1)
        #     @printf("Turn off generator %d . . .\n", neg_gen)
        # end
        # flush(stdout)

        opt = Option()
        opt.obj_gencost = true
        opt.has_ramping = true
        opt.freq_ctrl = true
        opt.sc_constr = true
        # opt.neg_g = [neg_gen]
        # opt.weight_freqctrl = weight_freqctrl

        demand = Load(load.pd[:,1:T], load.qd[:,1:T])
        m_sc = get_mpcmodel(circuit, demand; opt = opt)
        init_x(m_sc, circuit, demand; opt = opt)
        em = ExaModel(m_sc)
        result = ipopt(em)
        # setsolver(m_sc, IpoptSolver(option_file_name="ipopt.opt"))
        # stat = solve(m_sc)
        # if stat != :Optimal
        #     println("Stat is not optimal: ", stat)
        #     return
        # end

        println("\nOmega and Frequencies:")
        for t in 1:T
            @printf("%.6e    %.16f\n",
                    value(m_sc[:omega][t]), 60 + 60*value(m_sc[:omega][t]))
        end

        pos_gen = findall(x -> x != neg_gen, collect(1:length(circuit.gen)))
        pg_old_x = zeros(length(pos_gen))
        # for (j,g) in enumerate(pos_gen)
        #    pg_old_x[j] = getvalue(m_sc[:Pg_f][g])
        # end

        # circuit = getcircuit(case, baseMVA, ramp_scale;
        #                      neg_line=neg_line, neg_gen=neg_gen)
        # save_circuit(circuit, pg_old_x, [neg_gen], "solution_circuit_"*profname)
    else
        # circuit, pg_old_x, neg_g = load_circuit("solution_circuit_"*profname)
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
    opt.obj_gencost = true
    opt.freq_ctrl = true
    opt.sc_constr = false
    opt.piecewise = piecewise_cost
    opt.powerflow_solve = powerflow_solve
    opt.neg_g = Int[]
    opt.weight_freqctrl = weight_freqctrl

    demand = Load(load.pd[:,1:T], load.qd[:,1:T])
    m_cur = get_mpcmodel(circuit, demand; opt=opt)
    em = ExaModel(m_cur)
    result = ipopt(em)
    print(result)
    # setsolver(m_cur, IpoptSolver(option_file_name="ipopt.opt"))

    if load_sol == 1
        stat = :Optimal
        load_solution(m_cur, "solution_"*profname)
        save_rampinfo(m_cur, T, num_gens, "rampinfo_"*profname, "w"; circuit=circuit)
    else
        if opt.freq_ctrl
            for g in 1:num_gens
                set_lower_bound(m_cur[:Pg_ref][1,g], pg_old_x[g])
                set_upper_bound(m_cur[:Pg_ref][1,g], pg_old_x[g])
                set_start_value.(m_cur[:Pg_ref][1,g], pg_old_x[g])
            end
        else
            for g in 1:num_gens
                set_lower_bound(m_cur[:Pg][1,g], pg_old_x[g])
                set_upper_bound(m_cur[:Pg][1,g], pg_old_x[g])
                set_start_value.(m_cur[:Pg][1,g], pg_old_x[g])
            end
        end

        # JuMP.build(m_cur)
        # set_warm_from_sc(m_cur, m_sc, circuit, T; neg_g=neg_gen)
        # in_cur = internalmodel(m_cur)
        # MathProgBase.setwarmstart!(in_cur, in_cur.inner.x)
        # MathProgBase.optimize!(in_cur)

        em = ExaModel(m_cur)
        result = ipopt(em)
        print(result)

        # stat = MathProgBase.status(in_cur)
        # if stat != :Optimal
        #     println("Stat is not optimal: ", stat)
        #     return
        # end

        # m_cur.objVal = MathProgBase.getobjval(in_cur)
        # m_cur.colVal = MathProgBase.getsolution(in_cur)

        # save_solution(m_cur, "solution_"*profname)
        # save_rampinfo(m_cur, T, num_gens, "rampinfo_"*profname, "w"; circuit=circuit)
    end

    # single_objval = get_mpcobjectivevalue(m_cur, circuit, 1)
    # total_single_objval += single_objval
    # total_objval += single_objval
    # print_statistics(m_cur, circuit, T, opt)

    # ---------------------------------------------------------------------
    # Warm-start option for Ipopt.
    # ---------------------------------------------------------------------

    # m_prev = nothing
    # for h in 2:H
    #     println("\nh = h", h, "\n")
    #     flush(stdout)

    #     m_prev = m_cur
    #     inner_prev = internalmodel(m_prev).inner

    #     # -----------------------------------------------------------------
    #     # Build a model for the current horizon: [h,h+T-1]
    #     # -----------------------------------------------------------------

    #     demand = Load(load.pd[:,h:h+T-1], load.qd[:,h:h+T-1])
    #     if perturb != 0
    #         Random.seed!(0)
    #         factor = rand((-5,5))*(rand(Float64,1)*perturb)
    #         demand.pd .= demand.pd .+ factor.*demand.pd
    #         demand.qd .= demand.qd .+ factor.*demand.qd
    #     end

    #     m_cur = get_mpcmodel(circuit, demand; opt=opt)

    #     # -----------------------------------------------------------------
    #     # Change the lower/upper bounds on the first variable based on
    #     # the solution of the previous model.
    #     # -----------------------------------------------------------------

    #     Pg_prev = Pg_cur = nothing

    #     if opt.freq_ctrl
    #         Pg_prev = m_prev[:Pg_ref]
    #         Pg_cur = m_cur[:Pg_ref]
    #     else
    #         Pg_prev = m_prev[:Pg]
    #         Pg_cur = m_cur[:Pg]
    #     end

    #     start = linearindex(Pg_prev[1])
    #     for g=1:num_gens
    #         setlowerbound(Pg_cur[1,g],
    #                       max(gen[g].Pmin, inner_prev.x[start+g-1] - gen[g].ramp_agc))
    #         setupperbound(Pg_cur[1,g],
    #                       min(gen[g].Pmax, inner_prev.x[start+g-1] + gen[g].ramp_agc))
    #     end

    #     if T >= 2 && warmstart
    #         if h == 2
    #             rm("rampinfo_p1_"*profname*".txt", force=true)
    #         end

    #         # -------------------------------------------------------------
    #         # For a warm-start,
    #         #
    #         #  - set options for Ipopt.
    #         #  - build a model.
    #         #  - copy primal/dual values and constraints values.
    #         # -------------------------------------------------------------

    #         optname = "ipopt.op"*string(optfile_num)
    #         setsolver(m_cur, IpoptSolver(option_file_name=optname))
    #         JuMP.build(m_cur)
    #         set_warm(m_cur, m_prev, circuit, demand;
    #                  warm_type=warm_type, profname=profname, opt=opt)

    #         if qp == 1
    #             d, m_qp = mpc_qp(m_cur; optfile=optname)
    #             copyto!(m_cur.colVal, 1, m_qp.colVal, 1, MathProgBase.numvar(m_qp))
    #             m_cur.objVal = get_mpcobjectivevalue(m_cur, circuit, T; opt=opt)
    #             constr_viol = get_mpcconstrviolation(m_cur, d)

    #             @printf("QP approximation objective..................: %18.16e\n", m_cur.objVal)
    #             @printf("QP approximation constraint violation.......: %18.16e\n", constr_viol)

    #             copyvars(internalmodel(m_cur).inner, 1,
    #                      internalmodel(m_qp).inner, 1,
    #                      MathProgBase.numvar(m_qp))
    #             copyconstrs(internalmodel(m_cur).inner, 1,
    #                         internalmodel(m_qp).inner, 1,
    #                         MathProgBase.numconstr(m_qp))
    #         else
    #             # -------------------------------------------------------------
    #             # Call optimize!() for a warm-start instead of solve() because
    #             # solve() will build a new internal model that ignores the
    #             # multipliers.
    #             #
    #             # Note that MathProgBase.optimize!() doesn't change the original
    #             # model's solution status, such as objVal, colVal, and redCosts.
    #             # If we need to access them using the original model, we should
    #             # manually assign them.
    #             # -------------------------------------------------------------

    #             in_cur = internalmodel(m_cur)   # IpoptMathProgModel is returned.
    #             MathProgBase.setwarmstart!(in_cur, in_cur.inner.x)
    #             MathProgBase.optimize!(in_cur)
    #             stat = MathProgBase.status(in_cur)

    #             if stat != :Infeasible && stat != :Unbounded
    #                 m_cur.objVal = MathProgBase.getobjval(in_cur)
    #                 m_cur.colVal = MathProgBase.getsolution(in_cur)
    #             end
    #         end
    #     else

    #         # -------------------------------------------------------------
    #         # Cold-start
    #         # -------------------------------------------------------------

    #         stat = solve_cold(m_cur, circuit, demand; powerflow_solve = powerflow_solve)
    #     end

    #     if stat != :Optimal && stat != :UserLimit
    #         println("Stat is not optimal: ", stat)
    #         return
    #     end

    #     save_rampinfo(m_cur, T, num_gens, "rampinfo_"*profname; circuit=circuit)

    #     single_objval = get_mpcobjectivevalue(m_cur, circuit, 1)
    #     total_single_objval += single_objval
    #     total_objval += single_objval

    #     print_statistics(m_cur, circuit, T, opt)
    # end

    # # Add the remaining cost.
    # total_objval += getobjectivevalue(m_cur) - single_objval

    @printf("Total objective value..............: %18.6e\n", total_objval)
    @printf("Total single objective value.......: %18.6e\n", total_single_objval)
    
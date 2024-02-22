using Ipopt, NLPModelsIpopt
using JuMP
using ExaModels, MadNLP, MadNLPGPU
using Random
using DelimitedFiles, Printf, CUDA
include("mpc_data.jl")
include("mpc_model.jl")
include("mpc_2.jl")

case  = "data/case1354pegase"
scen = "data/case1354pegase/halfhour_30"
T = 5
H = 1
load_scale = 1.0
ramp_scale = 0.7
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
        # set_optimizer(m_cur, Ipopt.Optimizer)
        # optimize!(m_cur)
        em = ExaModel(m_cur, backend = CUDABackend())
        result = madnlp(em)
        #result = ipopt(em)
        print(result)
        

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
        # set_optimizer(m_cur, Ipopt.Optimizer)
        # optimize!(m_cur)
        em = ExaModel(m_sc, backend = CUDABackend())
        result = madnlp(em)
        #result = ipopt(em)
        # setsolver(m_sc, IpoptSolver(option_file_name="ipopt.opt"))
        # stat = solve(m_sc)
        # if stat != :Optimal
        #     println("Stat is not optimal: ", stat)
        #     return
        # end

        println("\nOmega and Frequencies:")
        # for t in 1:T
        #     @printf("%.6e    %.16f\n",
        #             value(m_sc[:omega][t]), 60 + 60*value(m_sc[:omega][t]))
        # end

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

   

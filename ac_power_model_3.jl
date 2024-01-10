function parse_ac_power_data(filename)
    data = PowerModels.parse_file(filename)
    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    arcdict = Dict(
        a=>k
        for (k,a) in enumerate(ref[:arcs]))
    busdict = Dict(
        k=>i
        for (i,(k,v)) in enumerate(ref[:bus]))
    gendict = Dict(
        k=>i
        for (i,(k,v)) in enumerate(ref[:gen]))
    branchdict = Dict(
        k=>i
        for (i,(k,v)) in enumerate(ref[:branch]))

    return (
        bus = [
            begin
                bus_loads = [ref[:load][l] for l in ref[:bus_loads][k]]
                bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][k]]
                pd = sum(load["pd"] for load in bus_loads; init = 0.)
                gs = sum(shunt["gs"] for shunt in bus_shunts; init = 0.)
                qd = sum(load["qd"] for load in bus_loads; init = 0.)
                bs = sum(shunt["bs"] for shunt in bus_shunts; init = 0.)
                (
                    i = busdict[k],
                    pd = pd, gs = gs, qd = qd, bs = bs
                )
            end
            for (k,v) in ref[:bus]],
        gen = [
            (i = gendict[k],
             cost1 = v["cost"][1], cost2 = v["cost"][2], cost3 = v["cost"][3], bus = busdict[v["gen_bus"]])
            for (k,v) in ref[:gen]],
        arc =[
            (i=k, rate_a = ref[:branch][l]["rate_a"], bus = busdict[i])
            for (k,(l,i,j)) in enumerate(ref[:arcs])],
        branch = [
            begin
                f_idx = arcdict[i, branch["f_bus"], branch["t_bus"]]
                t_idx = arcdict[i, branch["t_bus"], branch["f_bus"]]
                g, b = PowerModels.calc_branch_y(branch)
                tr, ti = PowerModels.calc_branch_t(branch)
                ttm = tr^2 + ti^2
                g_fr = branch["g_fr"]
                b_fr = branch["b_fr"]
                g_to = branch["g_to"]
                b_to = branch["b_to"]
                c1 = (-g*tr-b*ti)/ttm
                c2 = (-b*tr+g*ti)/ttm
                c3 = (-g*tr+b*ti)/ttm
                c4 = (-b*tr-g*ti)/ttm
                c5 = (g+g_fr)/ttm
                c6 = (b+b_fr)/ttm
                c7 = (g+g_to)
                c8 = (b+b_to)
                (
                    i = branchdict[i],
                    j = 1,
                    f_idx = f_idx,
                    t_idx = t_idx,
                    f_bus = busdict[branch["f_bus"]],
                    t_bus = busdict[branch["t_bus"]],
                    c1 = c1,
                    c2 = c2,
                    c3 = c3,
                    c4 = c4,
                    c5 = c5,
                    c6 = c6,
                    c7 = c7,
                    c8 = c8,
                    rate_a_sq = branch["rate_a"]^2,
                )
            end
            for (i,branch) in ref[:branch]],
        ref_buses = [busdict[i] for (i,k) in ref[:ref_buses]],
        vmax = [
            v["vmax"] for (k,v) in ref[:bus]],
        vmin = [
            v["vmin"] for (k,v) in ref[:bus]],
        pmax = [
            v["pmax"] for (k,v) in ref[:gen]],
        pmin = [
            v["pmin"] for (k,v) in ref[:gen]],
        qmax = [
            v["qmax"] for (k,v) in ref[:gen]],
        qmin = [
            v["qmin"] for (k,v) in ref[:gen]],
        rate_a =[
            ref[:branch][l]["rate_a"]
            for (k,(l,i,j)) in enumerate(ref[:arcs])],
        angmax = [b["angmax"] for (i,b) in ref[:branch]],
        angmin = [b["angmin"] for (i,b) in ref[:branch]],
    )
end

convert_data(data::N, backend) where {names, N <: NamedTuple{names}} = NamedTuple{names}(ExaModels.convert_array(d,backend) for d in data)

parse_ac_power_data(filename, backend) = convert_data(parse_ac_power_data(filename), backend)

function ac_power_model(
    filename;
    backend = nothing,
    T = Float64
    )

    data = parse_ac_power_data(filename, backend)

    w = ExaCore(T, backend)

 

    
    # w = ExaCore(T, backend)

    va = variable(
        w, 
        length(data.bus);
    )

    va_1 = variable(
        w, 
        length(data.bus);
    )

    va_2 = variable(
        w, 
        length(data.bus);
    )

    va_3 = variable(
        w, 
        length(data.bus);
    )

    va_4 = variable(
        w, 
        length(data.bus);
    )

    va_5 = variable(
        w, 
        length(data.bus);
    )

    va_6 = variable(
        w, 
        length(data.bus);
    )

    # va_7 = variable(
    #     w, 
    #     length(data.bus);
    # )

    # va_8 = variable(
    #     w, 
    #     length(data.bus);
    # )

    # va_9 = variable(
    #     w, 
    #     length(data.bus);
    # )

    
    vm = variable(
        w,
        length(data.bus);
        start = fill!(similar(data.bus,Float64),1.),
        lvar = data.vmin,
        uvar = data.vmax
    )


    

    vm_1 = variable(
        w,
        length(data.bus);
        start = fill!(similar(data.bus,Float64),1.),
        lvar = data.vmin,
        uvar = data.vmax
    )

    vm_2 = variable(
        w,
        length(data.bus);
        start = fill!(similar(data.bus,Float64),1.),
        lvar = data.vmin,
        uvar = data.vmax
    )

    vm_3 = variable(
        w,
        length(data.bus);
        start = fill!(similar(data.bus,Float64),1.),
        lvar = data.vmin,
        uvar = data.vmax
    )

    vm_4 = variable(
        w,
        length(data.bus);
        start = fill!(similar(data.bus,Float64),1.),
        lvar = data.vmin,
        uvar = data.vmax
    )


    vm_5 = variable(
        w,
        length(data.bus);
        start = fill!(similar(data.bus,Float64),1.),
        lvar = data.vmin,
        uvar = data.vmax
    )

    vm_6 = variable(
        w,
        length(data.bus);
        start = fill!(similar(data.bus,Float64),1.),
        lvar = data.vmin,
        uvar = data.vmax
    )

    # vm_7 = variable(
    #     w,
    #     length(data.bus);
    #     start = fill!(similar(data.bus,Float64),1.),
    #     lvar = data.vmin,
    #     uvar = data.vmax
    # )

    # vm_8 = variable(
    #     w,
    #     length(data.bus);
    #     start = fill!(similar(data.bus,Float64),1.),
    #     lvar = data.vmin,
    #     uvar = data.vmax
    # )

    # vm_9 = variable(
    #     w,
    #     length(data.bus);
    #     start = fill!(similar(data.bus,Float64),1.),
    #     lvar = data.vmin,
    #     uvar = data.vmax
    # )


    

    pg = variable(
        w,
        length(data.gen);
        lvar = data.pmin,
        uvar = data.pmax
    )



    pg_1 = variable(
        w,
        length(data.gen);
        lvar = data.pmin,
        uvar = data.pmax
    )

    pg_2 = variable(
        w,
        length(data.gen);
        lvar = data.pmin,
        uvar = data.pmax
    )

    pg_3 = variable(
        w,
        length(data.gen);
        lvar = data.pmin,
        uvar = data.pmax
    )

    pg_4 = variable(
        w,
        length(data.gen);
        lvar = data.pmin,
        uvar = data.pmax
    )

    pg_5 = variable(
        w,
        length(data.gen);
        lvar = data.pmin,
        uvar = data.pmax
    )

    pg_6 = variable(
        w,
        length(data.gen);
        lvar = data.pmin,
        uvar = data.pmax
    )

    # pg_7 = variable(
    #     w,
    #     length(data.gen);
    #     lvar = data.pmin,
    #     uvar = data.pmax
    # )

    # pg_8 = variable(
    #     w,
    #     length(data.gen);
    #     lvar = data.pmin,
    #     uvar = data.pmax
    # )

    # pg_9 = variable(
    #     w,
    #     length(data.gen);
    #     lvar = data.pmin,
    #     uvar = data.pmax
    # )

    qg = variable(
        w,
        length(data.gen);
        lvar = data.qmin,
        uvar = data.qmax
    )

    
    qg_1 = variable(
        w,
        length(data.gen);
        lvar = data.qmin,
        uvar = data.qmax
    )

    qg_2 = variable(
        w,
        length(data.gen);
        lvar = data.qmin,
        uvar = data.qmax
    )

    qg_3 = variable(
        w,
        length(data.gen);
        lvar = data.qmin,
        uvar = data.qmax
    )

    qg_4 = variable(
        w,
        length(data.gen);
        lvar = data.qmin,
        uvar = data.qmax
    )

    qg_5 = variable(
        w,
        length(data.gen);
        lvar = data.qmin,
        uvar = data.qmax
    )

    qg_6 = variable(
        w,
        length(data.gen);
        lvar = data.qmin,
        uvar = data.qmax
    )

    # qg_7 = variable(
    #     w,
    #     length(data.gen);
    #     lvar = data.qmin,
    #     uvar = data.qmax
    # )

    # qg_8 = variable(
    #     w,
    #     length(data.gen);
    #     lvar = data.qmin,
    #     uvar = data.qmax
    # )

    # qg_9 = variable(
    #     w,
    #     length(data.gen);
    #     lvar = data.qmin,
    #     uvar = data.qmax
    # )

   

    p = variable(
        w,
        length(data.arc);
        lvar = -data.rate_a,
        uvar = data.rate_a
    )

    
    p_1 = variable(
        w,
        length(data.arc);
        lvar = -data.rate_a,
        uvar = data.rate_a
    )

    p_2 = variable(
        w,
        length(data.arc);
        lvar = -data.rate_a,
        uvar = data.rate_a
    )

    p_3 = variable(
        w,
        length(data.arc);
        lvar = -data.rate_a,
        uvar = data.rate_a
    )

    p_4 = variable(
        w,
        length(data.arc);
        lvar = -data.rate_a,
        uvar = data.rate_a
    )

    p_5 = variable(
        w,
        length(data.arc);
        lvar = -data.rate_a,
        uvar = data.rate_a
    )

    p_6 = variable(
        w,
        length(data.arc);
        lvar = -data.rate_a,
        uvar = data.rate_a
    )

    # p_7 = variable(
    #     w,
    #     length(data.arc);
    #     lvar = -data.rate_a,
    #     uvar = data.rate_a
    # )

    # p_8 = variable(
    #     w,
    #     length(data.arc);
    #     lvar = -data.rate_a,
    #     uvar = data.rate_a
    # )

    # p_9 = variable(
    #     w,
    #     length(data.arc);
    #     lvar = -data.rate_a,
    #     uvar = data.rate_a
    # )


    q = variable(
        w,
        length(data.arc);
        lvar = -data.rate_a,
        uvar = data.rate_a
    )

    q_1 = variable(
        w,
        length(data.arc);
        lvar = -data.rate_a,
        uvar = data.rate_a
    )

    q_2 = variable(
        w,
        length(data.arc);
        lvar = -data.rate_a,
        uvar = data.rate_a
    )

    q_3 = variable(
        w,
        length(data.arc);
        lvar = -data.rate_a,
        uvar = data.rate_a
    )

    q_4 = variable(
        w,
        length(data.arc);
        lvar = -data.rate_a,
        uvar = data.rate_a
    )

    q_9 = variable(
        w,
        length(data.arc);
        lvar = -data.rate_a,
        uvar = data.rate_a
    )

    q_5 = variable(
        w,
        length(data.arc);
        lvar = -data.rate_a,
        uvar = data.rate_a
    )

    q_6 = variable(
        w,
        length(data.arc);
        lvar = -data.rate_a,
        uvar = data.rate_a
    )

    # q_7 = variable(
    #     w,
    #     length(data.arc);
    #     lvar = -data.rate_a,
    #     uvar = data.rate_a
    # )

    # q_8 = variable(
    #     w,
    #     length(data.arc);
    #     lvar = -data.rate_a,
    #     uvar = data.rate_a
    # )

    # q_9 = variable(
    #     w,
    #     length(data.arc);
    #     lvar = -data.rate_a,
    #     uvar = data.rate_a
    # )

    

    o = objective(
        w,
        g.cost1 * pg[g.i]^2 + g.cost2 * pg[g.i] + g.cost3
        for g in data.gen)

    c1 = constraint(
        w,
        va[i] for i in data.ref_buses)

   
    
    c1_1 = constraint(
        w,
        va_1[i] for i in data.ref_buses)

    c1_2 = constraint(
        w,
        va_2[i] for i in data.ref_buses)

    c1_3 = constraint(
        w,
        va_3[i] for i in data.ref_buses)
    
    c1_4 = constraint(
        w,
        va_4[i] for i in data.ref_buses)

    c1_5 = constraint(
        w,
        va_5[i] for i in data.ref_buses)
    
    c1_6 = constraint(
        w,
        va_6[i] for i in data.ref_buses)
    
    # c1_7 = constraint(
    #     w,
    #     va_7[i] for i in data.ref_buses)
        
    # c1_8 = constraint(
    #     w,
    #     va_8[i] for i in data.ref_buses)
    
    # c1_9 = constraint(
    #     w,
    #     va_9[i] for i in data.ref_buses)

    c2 = constraint(
        w,
        p[b.f_idx]
        - b.c5*vm[b.f_bus]^2
        - b.c3*(vm[b.f_bus]*vm[b.t_bus]*cos(va[b.f_bus]-va[b.t_bus]))
        - b.c4*(vm[b.f_bus]*vm[b.t_bus]*sin(va[b.f_bus]-va[b.t_bus]))
        for b in data.branch)
    

    c2_1 = constraint(
        w,
        p_1[b.f_idx]
        - b.c5*vm_1[b.f_bus]^2
        - b.c3*(vm_1[b.f_bus]*vm_1[b.t_bus]*cos(va_1[b.f_bus]-va_1[b.t_bus]))
        - b.c4*(vm_1[b.f_bus]*vm_1[b.t_bus]*sin(va_1[b.f_bus]-va_1[b.t_bus]))
        for b in data.branch)
    
    c2_2 = constraint(
        w,
        p_2[b.f_idx]
        - b.c5*vm_2[b.f_bus]^2
        - b.c3*(vm_2[b.f_bus]*vm_2[b.t_bus]*cos(va_2[b.f_bus]-va_2[b.t_bus]))
        - b.c4*(vm_2[b.f_bus]*vm_2[b.t_bus]*sin(va_2[b.f_bus]-va_2[b.t_bus]))
        for b in data.branch)
    
    c2_3 = constraint(
        w,
        p_3[b.f_idx]
        - b.c5*vm_3[b.f_bus]^2
        - b.c3*(vm_3[b.f_bus]*vm_3[b.t_bus]*cos(va_3[b.f_bus]-va_3[b.t_bus]))
        - b.c4*(vm_3[b.f_bus]*vm_3[b.t_bus]*sin(va_3[b.f_bus]-va_3[b.t_bus]))
        for b in data.branch)
    
    c2_4 = constraint(
        w,
        p_4[b.f_idx]
        - b.c5*vm_4[b.f_bus]^2
        - b.c3*(vm_4[b.f_bus]*vm_4[b.t_bus]*cos(va_4[b.f_bus]-va_4[b.t_bus]))
        - b.c4*(vm_4[b.f_bus]*vm_3[b.t_bus]*sin(va_4[b.f_bus]-va_4[b.t_bus]))
        for b in data.branch)

    c2_5 = constraint(
        w,
        p_5[b.f_idx]
        - b.c5*vm_5[b.f_bus]^2
        - b.c3*(vm_5[b.f_bus]*vm_5[b.t_bus]*cos(va_5[b.f_bus]-va_5[b.t_bus]))
        - b.c4*(vm_5[b.f_bus]*vm_5[b.t_bus]*sin(va_5[b.f_bus]-va_5[b.t_bus]))
        for b in data.branch)
    
    c2_6 = constraint(
        w,
        p_6[b.f_idx]
        - b.c5*vm_6[b.f_bus]^2
        - b.c3*(vm_6[b.f_bus]*vm_6[b.t_bus]*cos(va_6[b.f_bus]-va_6[b.t_bus]))
        - b.c4*(vm_6[b.f_bus]*vm_6[b.t_bus]*sin(va_6[b.f_bus]-va_6[b.t_bus]))
        for b in data.branch)
    
    # c2_7 = constraint(
    #     w,
    #     p_7[b.f_idx]
    #     - b.c5*vm_7[b.f_bus]^2
    #     - b.c3*(vm_7[b.f_bus]*vm_7[b.t_bus]*cos(va_7[b.f_bus]-va_7[b.t_bus]))
    #     - b.c4*(vm_7[b.f_bus]*vm_7[b.t_bus]*sin(va_7[b.f_bus]-va_7[b.t_bus]))
    #     for b in data.branch)
    
    # c2_8 = constraint(
    #     w,
    #     p_8[b.f_idx]
    #     - b.c5*vm_8[b.f_bus]^2
    #     - b.c3*(vm_8[b.f_bus]*vm_8[b.t_bus]*cos(va_8[b.f_bus]-va_8[b.t_bus]))
    #     - b.c4*(vm_8[b.f_bus]*vm_8[b.t_bus]*sin(va_8[b.f_bus]-va_8[b.t_bus]))
    #     for b in data.branch)    
    
    # c2_9 = constraint(
    #     w,
    #     p_9[b.f_idx]
    #     - b.c5*vm_9[b.f_bus]^2
    #     - b.c3*(vm_9[b.f_bus]*vm_9[b.t_bus]*cos(va_9[b.f_bus]-va_9[b.t_bus]))
    #     - b.c4*(vm_9[b.f_bus]*vm_9[b.t_bus]*sin(va_9[b.f_bus]-va_9[b.t_bus]))
    #     for b in data.branch)  
    

    c3 = constraint(
        w,
        q[b.f_idx]
        + b.c6*vm[b.f_bus]^2
        + b.c4*(vm[b.f_bus]*vm[b.t_bus]*cos(va[b.f_bus]-va[b.t_bus]))
        - b.c3*(vm[b.f_bus]*vm_1[b.t_bus]*sin(va[b.f_bus]-va[b.t_bus]))
        for b in data.branch)
    
    c3_1 = constraint(
        w,
        q_1[b.f_idx]
        + b.c6*vm_1[b.f_bus]^2
        + b.c4*(vm_1[b.f_bus]*vm_1[b.t_bus]*cos(va_1[b.f_bus]-va_1[b.t_bus]))
        - b.c3*(vm_1[b.f_bus]*vm_1[b.t_bus]*sin(va_1[b.f_bus]-va_1[b.t_bus]))
        for b in data.branch)
    
    c3_2 = constraint(
        w,
        q_2[b.f_idx]
        + b.c6*vm_2[b.f_bus]^2
        + b.c4*(vm_2[b.f_bus]*vm_2[b.t_bus]*cos(va_2[b.f_bus]-va_2[b.t_bus]))
        - b.c3*(vm_2[b.f_bus]*vm_2[b.t_bus]*sin(va_2[b.f_bus]-va_2[b.t_bus]))
        for b in data.branch)
    
    c3_3 = constraint(
        w,
        q_3[b.f_idx]
        + b.c6*vm_3[b.f_bus]^2
        + b.c4*(vm_3[b.f_bus]*vm_3[b.t_bus]*cos(va_3[b.f_bus]-va_3[b.t_bus]))
        - b.c3*(vm_3[b.f_bus]*vm_3[b.t_bus]*sin(va_3[b.f_bus]-va_3[b.t_bus]))
        for b in data.branch)

    c3_4 = constraint(
        w,
        q_4[b.f_idx]
        + b.c6*vm_4[b.f_bus]^2
        + b.c4*(vm_4[b.f_bus]*vm_4[b.t_bus]*cos(va_4[b.f_bus]-va_4[b.t_bus]))
        - b.c3*(vm_4[b.f_bus]*vm_4[b.t_bus]*sin(va_4[b.f_bus]-va_4[b.t_bus]))
        for b in data.branch)

    c3_5 = constraint(
        w,
        q_5[b.f_idx]
        + b.c6*vm_5[b.f_bus]^2
        + b.c4*(vm_5[b.f_bus]*vm_5[b.t_bus]*cos(va_5[b.f_bus]-va_5[b.t_bus]))
        - b.c3*(vm_5[b.f_bus]*vm_5[b.t_bus]*sin(va_5[b.f_bus]-va_5[b.t_bus]))
        for b in data.branch)
    
    c3_6 = constraint(
        w,
        q_6[b.f_idx]
        + b.c6*vm_6[b.f_bus]^2
        + b.c4*(vm_6[b.f_bus]*vm_6[b.t_bus]*cos(va_6[b.f_bus]-va_6[b.t_bus]))
        - b.c3*(vm_6[b.f_bus]*vm_6[b.t_bus]*sin(va_6[b.f_bus]-va_6[b.t_bus]))
        for b in data.branch)
    
    # c3_7 = constraint(
    #     w,
    #     q_7[b.f_idx]
    #     + b.c6*vm_7[b.f_bus]^2
    #     + b.c4*(vm_7[b.f_bus]*vm_7[b.t_bus]*cos(va_7[b.f_bus]-va_7[b.t_bus]))
    #     - b.c3*(vm_7[b.f_bus]*vm_7[b.t_bus]*sin(va_7[b.f_bus]-va_7[b.t_bus]))
    #     for b in data.branch)

    # c3_8 = constraint(
    #     w,
    #     q_8[b.f_idx]
    #     + b.c6*vm_8[b.f_bus]^2
    #     + b.c4*(vm_8[b.f_bus]*vm_8[b.t_bus]*cos(va_8[b.f_bus]-va_8[b.t_bus]))
    #     - b.c3*(vm_8[b.f_bus]*vm_8[b.t_bus]*sin(va_8[b.f_bus]-va_8[b.t_bus]))
    #     for b in data.branch)

    # c3_9 = constraint(
    #     w,
    #     q_9[b.f_idx]
    #     + b.c6*vm_9[b.f_bus]^2
    #     + b.c4*(vm_9[b.f_bus]*vm_9[b.t_bus]*cos(va_9[b.f_bus]-va_9[b.t_bus]))
    #     - b.c3*(vm_9[b.f_bus]*vm_9[b.t_bus]*sin(va_9[b.f_bus]-va_9[b.t_bus]))
    #     for b in data.branch)            
        
        
    c4 = constraint(
         w,
         p[b.t_idx]
         - b.c7*vm[b.t_bus]^2
         - b.c1*(vm[b.t_bus]*vm[b.f_bus]*cos(va[b.t_bus]-va[b.f_bus]))
         - b.c2*(vm[b.t_bus]*vm[b.f_bus]*sin(va[b.t_bus]-va[b.f_bus]))
         for b in data.branch
     )
    
    

    
    c4_1 = constraint(
        w,
        p_1[b.t_idx]
        - b.c7*vm_1[b.t_bus]^2
        - b.c1*(vm_1[b.t_bus]*vm_1[b.f_bus]*cos(va_1[b.t_bus]-va_1[b.f_bus]))
        - b.c2*(vm_1[b.t_bus]*vm_1[b.f_bus]*sin(va_1[b.t_bus]-va_1[b.f_bus]))
        for b in data.branch
    )


    c4_2 = constraint(
        w,
        p_2[b.t_idx]
        - b.c7*vm_2[b.t_bus]^2
        - b.c1*(vm_2[b.t_bus]*vm_2[b.f_bus]*cos(va_2[b.t_bus]-va_2[b.f_bus]))
        - b.c2*(vm_2[b.t_bus]*vm_2[b.f_bus]*sin(va_2[b.t_bus]-va_2[b.f_bus]))
        for b in data.branch
    )


    c4_3 = constraint(
        w,
        p_3[b.t_idx]
        - b.c7*vm_3[b.t_bus]^2
        - b.c1*(vm_3[b.t_bus]*vm_3[b.f_bus]*cos(va_3[b.t_bus]-va_3[b.f_bus]))
        - b.c2*(vm_3[b.t_bus]*vm_3[b.f_bus]*sin(va_3[b.t_bus]-va_3[b.f_bus]))
        for b in data.branch
    )


    c4_4 = constraint(
        w,
        p_2[b.t_idx]
        - b.c7*vm_4[b.t_bus]^2
        - b.c1*(vm_4[b.t_bus]*vm_4[b.f_bus]*cos(va_4[b.t_bus]-va_4[b.f_bus]))
        - b.c2*(vm_4[b.t_bus]*vm_4[b.f_bus]*sin(va_4[b.t_bus]-va_4[b.f_bus]))
        for b in data.branch
    )

    c4_5 = constraint(
        w,
        p_5[b.t_idx]
        - b.c7*vm_5[b.t_bus]^2
        - b.c1*(vm_5[b.t_bus]*vm_5[b.f_bus]*cos(va_5[b.t_bus]-va_5[b.f_bus]))
        - b.c2*(vm_5[b.t_bus]*vm_5[b.f_bus]*sin(va_5[b.t_bus]-va_5[b.f_bus]))
        for b in data.branch
    )


    c4_6 = constraint(
        w,
        p_6[b.t_idx]
        - b.c7*vm_6[b.t_bus]^2
        - b.c1*(vm_6[b.t_bus]*vm_6[b.f_bus]*cos(va_6[b.t_bus]-va_6[b.f_bus]))
        - b.c2*(vm_6[b.t_bus]*vm_6[b.f_bus]*sin(va_6[b.t_bus]-va_6[b.f_bus]))
        for b in data.branch
    )


    # c4_7 = constraint(
    #     w,
    #     p_7[b.t_idx]
    #     - b.c7*vm_7[b.t_bus]^2
    #     - b.c1*(vm_7[b.t_bus]*vm_7[b.f_bus]*cos(va_7[b.t_bus]-va_7[b.f_bus]))
    #     - b.c2*(vm_7[b.t_bus]*vm_7[b.f_bus]*sin(va_7[b.t_bus]-va_7[b.f_bus]))
    #     for b in data.branch
    # )


    # c4_8 = constraint(
    #     w,
    #     p_8[b.t_idx]
    #     - b.c7*vm_8[b.t_bus]^2
    #     - b.c1*(vm_8[b.t_bus]*vm_8[b.f_bus]*cos(va_8[b.t_bus]-va_8[b.f_bus]))
    #     - b.c2*(vm_8[b.t_bus]*vm_8[b.f_bus]*sin(va_8[b.t_bus]-va_8[b.f_bus]))
    #     for b in data.branch
    # )

    # c4_9 = constraint(
    #     w,
    #     p_9[b.t_idx]
    #     - b.c7*vm_9[b.t_bus]^2
    #     - b.c1*(vm_9[b.t_bus]*vm_9[b.f_bus]*cos(va_9[b.t_bus]-va_9[b.f_bus]))
    #     - b.c2*(vm_9[b.t_bus]*vm_9[b.f_bus]*sin(va_9[b.t_bus]-va_9[b.f_bus]))
    #     for b in data.branch
    # )  

    c5 = constraint(
         w,
         q[b.t_idx]
         + b.c8*vm[b.t_bus]^2
         + b.c2*(vm[b.t_bus]*vm[b.f_bus]*cos(va[b.t_bus]-va[b.f_bus]))
         - b.c1*(vm[b.t_bus]*vm[b.f_bus]*sin(va[b.t_bus]-va[b.f_bus]))
         for b in data.branch)

    c5_1 = constraint(
        w,
        q_1[b.t_idx]
        + b.c8*vm_1[b.t_bus]^2
        + b.c2*(vm_1[b.t_bus]*vm_1[b.f_bus]*cos(va_1[b.t_bus]-va_1[b.f_bus]))
        - b.c1*(vm_1[b.t_bus]*vm_1[b.f_bus]*sin(va_1[b.t_bus]-va_1[b.f_bus]))
        for b in data.branch)

    c5_2 = constraint(
        w,
        q_2[b.t_idx]
        + b.c8*vm_2[b.t_bus]^2
        + b.c2*(vm_2[b.t_bus]*vm_2[b.f_bus]*cos(va_2[b.t_bus]-va_2[b.f_bus]))
        - b.c1*(vm_2[b.t_bus]*vm_2[b.f_bus]*sin(va_2[b.t_bus]-va_2[b.f_bus]))
        for b in data.branch)

    c5_3 = constraint(
        w,
        q_3[b.t_idx]
        + b.c8*vm_3[b.t_bus]^2
        + b.c2*(vm_3[b.t_bus]*vm_3[b.f_bus]*cos(va_3[b.t_bus]-va_3[b.f_bus]))
        - b.c1*(vm_3[b.t_bus]*vm_3[b.f_bus]*sin(va_3[b.t_bus]-va_3[b.f_bus]))
        for b in data.branch)

    c5_4 = constraint(
        w,
        q_4[b.t_idx]
        + b.c8*vm_4[b.t_bus]^2
        + b.c2*(vm_4[b.t_bus]*vm_4[b.f_bus]*cos(va_4[b.t_bus]-va_4[b.f_bus]))
        - b.c1*(vm_4[b.t_bus]*vm_4[b.f_bus]*sin(va_4[b.t_bus]-va_4[b.f_bus]))
        for b in data.branch)

    c5_5 = constraint(
        w,
        q_5[b.t_idx]
        + b.c8*vm_5[b.t_bus]^2
        + b.c2*(vm_5[b.t_bus]*vm_5[b.f_bus]*cos(va_5[b.t_bus]-va_5[b.f_bus]))
        - b.c1*(vm_5[b.t_bus]*vm_5[b.f_bus]*sin(va_5[b.t_bus]-va_5[b.f_bus]))
        for b in data.branch)

    c5_6 = constraint(
        w,
        q_6[b.t_idx]
        + b.c8*vm_6[b.t_bus]^2
        + b.c2*(vm_6[b.t_bus]*vm_6[b.f_bus]*cos(va_6[b.t_bus]-va_6[b.f_bus]))
        - b.c1*(vm_6[b.t_bus]*vm_6[b.f_bus]*sin(va_6[b.t_bus]-va_6[b.f_bus]))
        for b in data.branch)

    # c5_7 = constraint(
    #     w,
    #     q_7[b.t_idx]
    #     + b.c8*vm_7[b.t_bus]^2
    #     + b.c2*(vm_7[b.t_bus]*vm_7[b.f_bus]*cos(va_7[b.t_bus]-va_7[b.f_bus]))
    #     - b.c1*(vm_7[b.t_bus]*vm_7[b.f_bus]*sin(va_7[b.t_bus]-va_7[b.f_bus]))
    #     for b in data.branch)

    # c5_8 = constraint(
    #     w,
    #     q_8[b.t_idx]
    #     + b.c8*vm_8[b.t_bus]^2
    #     + b.c2*(vm_8[b.t_bus]*vm_8[b.f_bus]*cos(va_8[b.t_bus]-va_8[b.f_bus]))
    #     - b.c1*(vm_8[b.t_bus]*vm_8[b.f_bus]*sin(va_8[b.t_bus]-va_8[b.f_bus]))
    #     for b in data.branch)

    # c5_9 = constraint(
    #     w,
    #     q_9[b.t_idx]
    #     + b.c8*vm_9[b.t_bus]^2
    #     + b.c2*(vm_9[b.t_bus]*vm_9[b.f_bus]*cos(va_9[b.t_bus]-va_9[b.f_bus]))
    #     - b.c1*(vm_9[b.t_bus]*vm_9[b.f_bus]*sin(va_9[b.t_bus]-va_9[b.f_bus]))
    #     for b in data.branch)             

    c6 = constraint(
        w,
        va[b.f_bus] - va[b.t_bus] for b in data.branch;
        lcon = data.angmin,
        ucon = data.angmax
    )
    

    c6_1 = constraint(
        w,
        va_1[b.f_bus] - va_1[b.t_bus] for b in data.branch;
            lcon = data.angmin,
            ucon = data.angmax
            )
    
    c6_2 = constraint(
        w,
        va_2[b.f_bus] - va_2[b.t_bus] for b in data.branch;
            lcon = data.angmin,
            ucon = data.angmax
            )


    c6_3 = constraint(
        w,
        va_3[b.f_bus] - va_3[b.t_bus] for b in data.branch;
            lcon = data.angmin,
            ucon = data.angmax
            )

    c6_4 = constraint(
        w,
        va_4[b.f_bus] - va_4[b.t_bus] for b in data.branch;
            lcon = data.angmin,
            ucon = data.angmax
            )

    c6_5 = constraint(
        w,
        va_5[b.f_bus] - va_5[b.t_bus] for b in data.branch;
        lcon = data.angmin,
        ucon = data.angmax
    )
    

    c6_6 = constraint(
        w,
        va_6[b.f_bus] - va_6[b.t_bus] for b in data.branch;
            lcon = data.angmin,
            ucon = data.angmax
            )
    
    # c6_7 = constraint(
    #     w,
    #     va_7[b.f_bus] - va_7[b.t_bus] for b in data.branch;
    #         lcon = data.angmin,
    #         ucon = data.angmax
    #         )


    # c6_8 = constraint(
    #     w,
    #     va_8[b.f_bus] - va_8[b.t_bus] for b in data.branch;
    #         lcon = data.angmin,
    #         ucon = data.angmax
    #         )

    # c6_9 = constraint(
    #     w,
    #     va_9[b.f_bus] - va_9[b.t_bus] for b in data.branch;
    #         lcon = data.angmin,
    #         ucon = data.angmax
    #         )
    
    
    c7 = constraint(
        w,
        p[b.f_idx]^2 + q[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c7_1 = constraint(
        w,
        p_1[b.f_idx]^2 + q_1[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c7_2 = constraint(
        w,
        p_2[b.f_idx]^2 + q_2[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c7_3 = constraint(
        w,
        p_3[b.f_idx]^2 + q_3[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c7_4 = constraint(
        w,
        p_4[b.f_idx]^2 + q_4[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )


    c7_5 = constraint(
        w,
        p_5[b.f_idx]^2 + q_5[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c7_6 = constraint(
        w,
        p_6[b.f_idx]^2 + q_6[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    # c7_7 = constraint(
    #     w,
    #     p_7[b.f_idx]^2 + q_7[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
    #     lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    # )

    # c7_8 = constraint(
    #     w,
    #     p_8[b.f_idx]^2 + q_8[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
    #     lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    # )

    # c7_9 = constraint(
    #     w,
    #     p_9[b.f_idx]^2 + q_9[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
    #     lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    # )
    
    c8 = constraint(
        w,
        p[b.t_idx]^2 + q[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c8_1 = constraint(
        w,
        p_1[b.t_idx]^2 + q_1[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c8_2 = constraint(
        w,
        p_2[b.t_idx]^2 + q_2[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c8_3 = constraint(
        w,
        p_3[b.t_idx]^2 + q_3[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c8_4 = constraint(
        w,
        p_4[b.t_idx]^2 + q_4[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c8_5 = constraint(
        w,
        p_5[b.t_idx]^2 + q_5[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c8_6 = constraint(
        w,
        p_6[b.t_idx]^2 + q_6[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    # c8_7 = constraint(
    #     w,
    #     p_7[b.t_idx]^2 + q_7[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
    #     lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    # )

    # c8_8 = constraint(
    #     w,
    #     p_8[b.t_idx]^2 + q_8[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
    #     lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    # )

    # c8_9 = constraint(
    #     w,
    #     p_9[b.t_idx]^2 + q_9[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
    #     lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    # )


    c9 = constraint(
        w,
        + b.pd
        + b.gs * vm[b.i]^2
        for b in data.bus)


    c9_1 = constraint(
        w,
        + b.pd
        + b.gs * vm_1[b.i]^2
        for b in data.bus)
    
    c9_2 = constraint(
        w,
        + b.pd
        + b.gs * vm_2[b.i]^2
        for b in data.bus)
    
    c9_3 = constraint(
        w,
        + b.pd
        + b.gs * vm_3[b.i]^2
        for b in data.bus)
    
    c9_4 = constraint(
        w,
        + b.pd
        + b.gs * vm_4[b.i]^2
        for b in data.bus)


    c9_5 = constraint(
        w,
        + b.pd
        + b.gs * vm_5[b.i]^2
        for b in data.bus)
    
    c9_6 = constraint(
        w,
        + b.pd
        + b.gs * vm_6[b.i]^2
        for b in data.bus)
    
    # c9_7 = constraint(
    #     w,
    #     + b.pd
    #     + b.gs * vm_7[b.i]^2
    #     for b in data.bus)
    
    # c9_8 = constraint(
    #     w,
    #     + b.pd
    #     + b.gs * vm_8[b.i]^2
    #     for b in data.bus)

    # c9_9 = constraint(
    #     w,
    #     + b.pd
    #     + b.gs * vm_9[b.i]^2
    #     for b in data.bus)

    c10 = constraint(
        w,
        + b.qd
        - b.bs * vm[b.i]^2
        for b in data.bus)
    
    c10_1 = constraint(
        w,
        + b.qd
        - b.bs * vm_1[b.i]^2
        for b in data.bus)


    c10_2 = constraint(
        w,
        + b.qd
        - b.bs * vm_2[b.i]^2
        for b in data.bus)
    
    c10_3 = constraint(
        w,
        + b.qd
        - b.bs * vm_3[b.i]^2
        for b in data.bus)
    
    c10_4 = constraint(
        w,
        + b.qd
        - b.bs * vm_4[b.i]^2
        for b in data.bus)


    c10_5 = constraint(
        w,
        + b.qd
        - b.bs * vm_5[b.i]^2
        for b in data.bus)
    
    c10_6 = constraint(
        w,
        + b.qd
        - b.bs * vm_6[b.i]^2
        for b in data.bus)


    # c10_7 = constraint(
    #     w,
    #     + b.qd
    #     - b.bs * vm_7[b.i]^2
    #     for b in data.bus)
    
    # c10_8 = constraint(
    #     w,
    #     + b.qd
    #     - b.bs * vm_8[b.i]^2
    #     for b in data.bus)
    
    # c10_9 = constraint(
    #     w,
    #     + b.qd
    #     - b.bs * vm_9[b.i]^2
    #     for b in data.bus)

    c11 = constraint!(
        w,
        c9,
        a.bus => p[a.i]
        for a in data.arc)
    
    c11_1 = constraint!(
        w,
        c9_1,
        a.bus => p_1[a.i]
        for a in data.arc)
    
    c11_2 = constraint!(
        w,
        c9_2,
        a.bus => p_2[a.i]
        for a in data.arc)
    
    c11_3 = constraint!(
        w,
        c9_3,
        a.bus => p_3[a.i]
        for a in data.arc)
    
    c11_4 = constraint!(
        w,
        c9_4,
        a.bus => p_4[a.i]
        for a in data.arc)

    c11_5 = constraint!(
        w,
        c9_5,
        a.bus => p_5[a.i]
        for a in data.arc)
    
    c11_6 = constraint!(
        w,
        c9_6,
        a.bus => p_6[a.i]
        for a in data.arc)
    
    # c11_7 = constraint!(
    #     w,
    #     c9_7,
    #     a.bus => p_7[a.i]
    #     for a in data.arc)
    
    # c11_8 = constraint!(
    #     w,
    #     c9_8,
    #     a.bus => p_8[a.i]
    #     for a in data.arc)
    
    # c11_9 = constraint!(
    #     w,
    #     c9_9,
    #     a.bus => p_9[a.i]
    #     for a in data.arc)


    c12 = constraint!(
        w,
        c10,
        a.bus => q[a.i]
        for a in data.arc)

    c12_1 = constraint!(
        w,
        c10_1,
        a.bus => q_1[a.i]
        for a in data.arc)

    
    c12_2 = constraint!(
        w,
        c10_2,
        a.bus => q_2[a.i]
        for a in data.arc)

    c12_3 = constraint!(
        w,
        c10_3,
        a.bus => q_3[a.i]
        for a in data.arc)
    
    c12_4 = constraint!(
        w,
        c10_4,
        a.bus => q_4[a.i]
        for a in data.arc)

    c12_5 = constraint!(
        w,
        c10_5,
        a.bus => q_6[a.i]
        for a in data.arc)

    c12_6 = constraint!(
        w,
        c10_6,
        a.bus => q_6[a.i]
        for a in data.arc)

    
    # c12_7 = constraint!(
    #     w,
    #     c10_7,
    #     a.bus => q_7[a.i]
    #     for a in data.arc)

    # c12_8 = constraint!(
    #     w,
    #     c10_3,
    #     a.bus => q_8[a.i]
    #     for a in data.arc)
    
    # c12_9 = constraint!(
    #     w,
    #     c10_4,
    #     a.bus => q_9[a.i]
    #     for a in data.arc)

    c13 = constraint!(
        w,
        c9,
        g.bus =>-pg[g.i]
        for g in data.gen)
        
    c13_1 = constraint!(
        w,
        c9_1,
        g.bus =>-pg_1[g.i]
        for g in data.gen)

    c13_2 = constraint!(
        w,
        c9_2,
        g.bus =>-pg_2[g.i]
        for g in data.gen)

    c13_3 = constraint!(
        w,
        c9_3,
        g.bus =>-pg_3[g.i]
        for g in data.gen)
    
    c13_4 = constraint!(
        w,
        c9_4,
        g.bus =>-pg_4[g.i]
        for g in data.gen)

    c13_5 = constraint!(
        w,
        c9_5,
        g.bus =>-pg_5[g.i]
        for g in data.gen)
        
    c13_6 = constraint!(
        w,
        c9_6,
        g.bus =>-pg_6[g.i]
        for g in data.gen)

    # c13_7 = constraint!(
    #     w,
    #     c9_7,
    #     g.bus =>-pg_7[g.i]
    #     for g in data.gen)

    # c13_8 = constraint!(
    #     w,
    #     c9_8,
    #     g.bus =>-pg_8[g.i]
    #     for g in data.gen)
    
    # c13_9 = constraint!(
    #     w,
    #     c9_9,
    #     g.bus =>-pg_9[g.i]
    #     for g in data.gen)


    c14 = constraint!(
        w,
        c10,
        g.bus =>-qg[g.i]
        for g in data.gen)
    
    c14_1 = constraint!(
        w,
        c10_1,
        g.bus =>-qg_1[g.i]
        for g in data.gen)
    
    c14_2 = constraint!(
        w,
        c10_2,
        g.bus =>-qg_2[g.i]
        for g in data.gen)

    c14_3 = constraint!(
        w,
        c10_3,
        g.bus =>-qg_3[g.i]
        for g in data.gen)
    
    c14_4 = constraint!(
        w,
        c10_4,
        g.bus =>-qg_4[g.i]
        for g in data.gen)
    
    c14_5 = constraint!(
        w,
        c10_5,
        g.bus =>-qg_5[g.i]
        for g in data.gen)
    
    c14_6 = constraint!(
        w,
        c10_6,
        g.bus =>-qg_6[g.i]
        for g in data.gen)
    
    # c14_7 = constraint!(
    #     w,
    #     c10_7,
    #     g.bus =>-qg_7[g.i]
    #     for g in data.gen)

    # c14_8 = constraint!(
    #     w,
    #     c10_8,
    #     g.bus =>-qg_8[g.i]
    #     for g in data.gen)
    
    # c14_9 = constraint!(
    #     w,
    #     c10_9,
    #     g.bus =>-qg_9[g.i]
    #     for g in data.gen)
    
    # c15_ramping = constraint(
    #     w,
    #      pg[g.i] - pg_1[g.i] <= 0.2* pg[g.i].uvar
    #     for g in data.gen)
    # c15_ramping_1 = constraint(
    #     w,
    #     pg[g.i] - pg_1[g.i] >= -0.2* pg[g.i].uvar
    #    for g in data.gen
    # ) 

    # c15_ramping = constraint(
    #     w,
    #     lcon = (-0.2*data.pmax[1, g.i] for g in data.gen),
    #     ucon = (0.2*data.pmax[1, g.i] for g in data.gen),
    #     pg[t,g.i] - pg[t-1, g.i] for g in data.gen;
    #     for (t,g) in itr0
    # )
    
    c15_ramping = constraint(
        w,
        pg_1[g.i] - pg[g.i] for g in data.gen;
        lcon = -0.2*data.pmax,
        ucon = 0.2*data.pmax
    )

    c15_ramping_1 = constraint(
        w,
        pg_2[g.i] - pg_1[g.i] for g in data.gen;
        lcon = -0.2*data.pmax,
        ucon = 0.2*data.pmax
    )

    c15_ramping_2 = constraint(
        w,
        pg_3[g.i] - pg_2[g.i] for g in data.gen;
        lcon = -0.2*data.pmax,
        ucon = 0.2*data.pmax
    )

    c15_ramping_3 = constraint(
        w,
        pg_4[g.i] - pg_3[g.i] for g in data.gen;
        lcon = -0.2*data.pmax,
        ucon = 0.2*data.pmax
    )


    c15_ramping_4 = constraint(
        w,
        pg_5[g.i] - pg_4[g.i] for g in data.gen;
        lcon = -0.2*data.pmax,
        ucon = 0.2*data.pmax
    )

    c15_ramping_5 = constraint(
        w,
        pg_6[g.i] - pg_5[g.i] for g in data.gen;
        lcon = -0.2*data.pmax,
        ucon = 0.2*data.pmax
    )

    # c15_ramping_6 = constraint(
    #     w,
    #     pg_7[g.i] - pg_6[g.i] for g in data.gen;
    #     lcon = -0.2*data.pmax,
    #     ucon = 0.2*data.pmax
    # )

    # c15_ramping_7 = constraint(
    #     w,
    #     pg_8[g.i] - pg_7[g.i] for g in data.gen;
    #     lcon = -0.2*data.pmax,
    #     ucon = 0.2*data.pmax
    # )

    # c15_ramping_8 = constraint(
    #     w,
    #     pg_9[g.i] - pg_8[g.i] for g in data.gen;
    #     lcon = -0.2*data.pmax,
    #     ucon = 0.2*data.pmax
    # )

    return ExaModel(w)

end
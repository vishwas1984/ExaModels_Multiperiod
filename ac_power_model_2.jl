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
        1:2,
        length(data.bus);
    )

    va_1 = variable(
        w, 
        1:2,
        length(data.bus);
    )

    
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

    return ExaModel(w)

end
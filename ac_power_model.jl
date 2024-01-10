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
    backend= nothing,
    T = Float64
    )

    data = parse_ac_power_data(filename, backend)

    w = ExaCore(T, backend)

 

    itr0 = ExaModels.convert_array(collect(Iterators.product(1:2, data.gen)), backend)
    itr1 = ExaModels.convert_array(collect(Iterators.product(2:2, data.gen)), backend)
    itr_vm1 = ExaModels.convert_array(collect(Iterators.product(1:2, data.branch)), backend)
    # itr_va1 = ExaModels.convert_array(collect(Iterators.product(1:2, data.branch)), backend)
    # itr_p1 = ExaModels.convert_array(collect(Iterators.product(1:2, data.branch)), backend)
    # itr_qg = ExaModels.convert_array(collect(Iterators.product(1:2, data.branch)), backend)
    itr_bus = ExaModels.convert_array(collect(Iterators.product(1:2, data.bus)), backend)
    itr_refbus = ExaModels.convert_array(collect(Iterators.product(1:2, data.ref_buses)), backend)
    itr_vmcon = ExaModels.convert_array(collect(Iterators.product(1:2, data.branch)), backend)
    # itr_vm2 = ExaModels.convert_array(collect(Iterators.product(2:2, data.bus)), backend)
    # itr1 = ExaModels.convert_array(collect(Iterators.product(1:2, 1:length(data.bus))), backend)

    # w = ExaCore(T, backend)

    va = variable(
        w, 
        1:2,
        length(data.bus);
    )

    # va_1 = variable(
    #     w, 
    #     1:2,
    #     length(data.bus);
    # )


    vm = variable(
        w,
        1:2,
        length(data.bus);
        start = fill!(similar(data.bus,Float64),1.),
        lvar = [p for i in 1:2, p in data.vmin],
        uvar = [p for i in 1:2, p in data.vmax]
    )

    # vm_1 = variable(
    #     w,
    #     length(data.bus);
    #     start = fill!(similar(data.bus,Float64),1.),
    #     lvar = data.vmin,
    #     uvar = data.vmax
    # )
    pg = variable(
        w,
        1:2,
        length(data.gen);
        lvar = [p for i in 1:2, p in data.pmin],
        uvar = [p for i in 1:2, p in data.pmax]
    )

    # pg_1 = variable(
    #     w,
    #     length(data.gen);
    #     lvar = data.pmin,
    #     uvar = data.pmax
    # )

    

    qg = variable(
        w,
        1:2,
        length(data.gen);
        lvar = [p for i in 1:2, p in data.qmin],
        uvar = [p for i in 1:2, p in data.qmax]
    )

    # qg_1 = variable(
    #     w,
    #     length(data.gen);
    #     lvar = data.qmin,
    #     uvar = data.qmax
    # )

    p = variable(
        w,
        1:2,
        length(data.arc);
        lvar = [p for i in 1:2, p in -data.rate_a],
        uvar = [p for i in 1:2, p in data.rate_a]
        # lvar = -data.rate_a,
        # uvar = data.rate_a
    )

    # p_1 = variable(
    #     w,
    #     length(data.arc);
    #     lvar = -data.rate_a,
    #     uvar = data.rate_a
    # )

    # q = variable(
    #     w,
    #     length(data.arc);
    #     lvar = -data.rate_a,
    #     uvar = data.rate_a
    # )


    q = variable(
        w,
        1:2,
        length(data.arc);
        lvar = [p for i in 1:2, p in -data.rate_a],
        uvar = [p for i in 1:2, p in data.rate_a]
        # lvar = -data.rate_a,
        # uvar = data.rate_a
    )



    # q_1 = variable(
    #     w,
    #     length(data.arc);
    #     lvar = -data.rate_a,
    #     uvar = data.rate_a
    # )

    o = objective(
        w,
        g.cost1 * pg[1,g.i]^2 + g.cost2 * pg[1,g.i] + g.cost3
        for g in data.gen)

    # c1 = constraint(
    #     w,
    #     va[i] for i in data.ref_buses)

    c1 = constraint(
        w,
        va[t,i] for (t,i) in itr_refbus
    )
    
    # c1_1 = constraint(
    #     w,
    #     va_1[i] for i in data.ref_buses)

    # c2 = constraint(
    #     w,
    #     p[b.f_idx]
    #     - b.c5*vm[b.f_bus]^2
    #     - b.c3*(vm[b.f_bus]*vm[b.t_bus]*cos(va[b.f_bus]-va[b.t_bus]))
    #     - b.c4*(vm[b.f_bus]*vm[b.t_bus]*sin(va[b.f_bus]-va[b.t_bus]))
    #     for b in data.branch)

    c2 = constraint(
        w,
        p[t, b.f_idx]
        - b.c5*vm[t,b.f_bus]^2
        - b.c3*(vm[t, b.f_bus]*vm[t, b.t_bus]*cos(va[t, b.f_bus]-va[t, b.t_bus]))
        - b.c4*(vm[t, b.f_bus]*vm[t, b.t_bus]*sin(va[t, b.f_bus]-va[t, b.t_bus]))
        for (t,b) in itr_vm1)



    # c2_1 = constraint(
    #     w,
    #     p_1[b.f_idx]
    #     - b.c5*vm_1[b.f_bus]^2
    #     - b.c3*(vm_1[b.f_bus]*vm_1[b.t_bus]*cos(va_1[b.f_bus]-va_1[b.t_bus]))
    #     - b.c4*(vm_1[b.f_bus]*vm_1[b.t_bus]*sin(va_1[b.f_bus]-va_1[b.t_bus]))
    #     for b in data.branch)

    
    
    

    c3 = constraint(
        w,
        q[t, b.f_idx]
        + b.c6*vm[t, b.f_bus]^2
        + b.c4*(vm[t, b.f_bus]*vm[t, b.t_bus]*cos(va[t, b.f_bus]-va[t, b.t_bus]))
        - b.c3*(vm[t, b.f_bus]*vm[t, b.t_bus]*sin(va[t, b.f_bus]-va[t, b.t_bus]))
        for (t,b) in itr_vm1)
    
    # c3_1 = constraint(
    #     w,
    #     q_1[b.f_idx]
    #     + b.c6*vm_1[b.f_bus]^2
    #     + b.c4*(vm_1[b.f_bus]*vm_1[b.t_bus]*cos(va_1[b.f_bus]-va_1[b.t_bus]))
    #     - b.c3*(vm_1[b.f_bus]*vm_1[b.t_bus]*sin(va_1[b.f_bus]-va_1[b.t_bus]))
    #     for b in data.branch)

    c4 = constraint(
        w,
        p[t, b.t_idx]
        - b.c7*vm[t, b.t_bus]^2
        - b.c1*(vm[t, b.t_bus]*vm[t, b.f_bus]*cos(va[t, b.t_bus]-va[t, b.f_bus]))
        - b.c2*(vm[t, b.t_bus]*vm[t, b.f_bus]*sin(va[t, b.t_bus]-va[t, b.f_bus]))
        for (t, b) in itr_vm1)
    
    # c4_1 = constraint(
    #     w,
    #     p_1[b.t_idx]
    #     - b.c7*vm_1[b.t_bus]^2
    #     - b.c1*(vm_1[b.t_bus]*vm_1[b.f_bus]*cos(va_1[b.t_bus]-va_1[b.f_bus]))
    #     - b.c2*(vm_1[b.t_bus]*vm_1[b.f_bus]*sin(va_1[b.t_bus]-va_1[b.f_bus]))
    #     for b in data.branch
    # )

    c5 = constraint(
        w,
        q[t, b.t_idx]
        + b.c8*vm[t, b.t_bus]^2
        + b.c2*(vm[t, b.t_bus]*vm[t, b.f_bus]*cos(va[t, b.t_bus]-va[t, b.f_bus]))
        - b.c1*(vm[t, b.t_bus]*vm[t, b.f_bus]*sin(va[t, b.t_bus]-va[t, b.f_bus]))
        for (t,b) in itr_vm1)

    # c5_1 = constraint(
    #     w,
    #     q_1[b.t_idx]
    #     + b.c8*vm_1[b.t_bus]^2
    #     + b.c2*(vm_1[b.t_bus]*vm_1[b.f_bus]*cos(va_1[b.t_bus]-va_1[b.f_bus]))
    #     - b.c1*(vm_1[b.t_bus]*vm_1[b.f_bus]*sin(va_1[b.t_bus]-va_1[b.f_bus]))
    #     for b in data.branch)
    


    c6 = constraint(
        w,
        va[t, b.f_bus] - va[t, b.t_bus] for (t,b) in itr_vm1;
            # lcon = data.angmin,
            # ucon = data.angmax
            lcon = [data.angmin[g.i] for (t,g) in itr_vm1],
            ucon = [data.angmax[g.i] for (t,g) in itr_vm1]
            )

    # c6_1 = constraint(
    #     w,
    #     va_1[b.f_bus] - va_1[b.t_bus] for b in data.branch;
    #         lcon = data.angmin,
    #         ucon = data.angmax
    #         )
    
    # c7 = constraint(
    #     w,
    #     p[1, b.f_idx]^2 + q[1, b.f_idx]^2 - b.rate_a_sq for (t,b) in itr_vm1;
    #     lcon = fill!(similar(itr_vmcon), -Inf)
    #          )


    # c7_lcon = vcat(
    #     (
    #         fill!(similar(data.branch, Float64, length(data.branch)), -Inf)
    #         for t = 1:2
    #     )...
    # )

    
    c7 = constraint(
        w,
        p[t, b.f_idx]^2 + q[t, b.f_idx]^2 - b.rate_a_sq for (t,b) in itr_vm1;
            lcon = -Inf
            )
    
    # c7 = constraint(
    #     w,
    #     p[2, b.f_idx]^2 + q[2, b.f_idx]^2 - b.rate_a_sq for b in data.branch;
    #         lcon = -Inf
    #         )
    
    # c8 = constraint(
    #     w,
    #     p[t, b.t_idx]^2 + q[t, b.t_idx]^2 - b.rate_a_sq for (t,b) in itr_vm1;
    #         # lcon = fill!(similar(data.branch,Float64,length(data.branch)),-Inf)
    #         )
    
    c8 = constraint(
        w,
        p[t, b.t_idx]^2 + q[t, b.t_idx]^2 - b.rate_a_sq for (t,b) in itr_vm1;
            lcon = -Inf
            )
    
    # c8_2 = constraint(
    #     w,
    #     p[2, b.t_idx]^2 + q[2, b.t_idx]^2 - b.rate_a_sq for b in data.branch;
    #     lcon = -Inf
    # )    

    c9 = constraint(
        w,
        + b.pd
        + b.gs * vm[t, b.i]^2
        for (t,b) in itr_bus)


    # c9_1 = constraint(
    #     w,
    #     + b.pd
    #     + b.gs * vm_1[b.i]^2
    #     for b in data.bus)

    c10 = constraint(
        w,
        + b.qd
        - b.bs * vm[t, b.i]^2
        for (t,b) in itr_bus)
    
    # c10_1 = constraint(
    #     w,
    #     + b.qd
    #     - b.bs * vm_1[b.i]^2
    #     for b in data.bus)


    itr_c11 = [(t,a.bus + (t-1)*length(data.bus), a.i) for t=1:2, a in data.arc]

    c11 = constraint!(
        w,
        c9,
        cidx => p[t,i]
        for (t,cidx,i) in itr_c11)
    
    # c11 = constraint!(
    #     w,
    #     c9,
    #     a.bus => p[t,a.i]
    #     for a in data.arc)
    
    # c11_1 = constraint!(
    #     w,
    #     c9_1,
    #     a.bus => p_1[t,a.i]
    #     for a in data.arc)

    c12 = constraint!(
        w,
        c10,
        cidx => q[t,i]
        for (t,cidx,i) in itr_c11)

    # c12 = constraint!(
    #     w,
    #     c10,
    #     a.bus => q[a.i]
    #     for a in data.arc)

    # c12_1 = constraint!(
    #     w,
    #     c10_1,
    #     a.bus => q_1[a.i]
    #     for a in data.arc)

    itr_c11 = [(t,g.bus + (t-1)*length(data.bus), g.i) for t=1:2, g in data.gen]

c13 = constraint!(
        w,
        c9,
        cidx =>-pg[t, i]
        for (t,cidx,i) in itr_c11)

# c13 = constraint!(
    #     w,
    #     c9,
    #     g.bus =>-pg[1, g.i]
    #     for g in data.gen)
        
    # c13_1 = constraint!(
    #     w,
    #     c9_1,
    #     g.bus =>-pg[2, g.i]
    #     for g in data.gen)

    c14 = constraint!(
        w,
        c10,
        cidx =>-qg[t,i]
        for (t,cidx,i) in itr_c11)
    # c14 = constraint!(
    #     w,
    #     c10,
    #     g.bus =>-qg[g.i]
    #     for g in data.gen)
    
    # c14_1 = constraint!(
    #     w,
    #     c10_1,
    #     g.bus =>-qg_1[g.i]
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
    
    # c15_ramping = constraint(
    #     w,
    #     pg[t,g.i] - pg[t-1, g.i] for (t,g) in itr1;
    #     lcon = [-0.2*data.pmax[g.i] for (t,g) in itr1],
    #     ucon = [0.2*data.pmax[g.i] for (t,g) in itr1],
    # )

    return ExaModel(w)

end

# using Pkg
using JuMP
using PowerModels
using PGLib
using Ipopt

# data = pglib("pglib_opf_case14_ieee.m")
# data = pglib("pglib_opf_case118_ieee.m")
# data = pglib("pglib_opf_case57_ieee.m")
# data = pglib("pglib_opf_case5_pjm.m")
data = pglib("pglib_opf_case4837_goc.m")

# Add zeros to turn linear objective functions into quadratic ones
# so that additional parameter checks are not required
PowerModels.standardize_cost_terms!(data, order=2)

# Adds reasonable rate_a values to branches without them
PowerModels.calc_thermal_limits!(data)

# use build_ref to filter out inactive components
data = PowerModels.build_ref(data)[:it][:pm][:nw][0]

function define_regions(network_data)
    # Set of all regions in network
    zones = Set([bus["zone"] for (i, bus) in data[:bus]])
    
    # Assigning each bus to a region
    buses = Dict{Int, Vector{Int64}}()
    for i in zones
        buses[i] = Vector{Int64}()
    end
    for (i, bus) in data[:bus]
        region = bus["zone"]
        push!(buses[region], i)
    end
    
    # Assigning each branch to a region
    branches = Dict{Int, Vector{Int64}}()
    coupled  = Vector{Int64}()

    for i in zones
        branches[i] = Vector{Int64}()
    end
    for (l, branch) in data[:branch]
        bus_fr = branch["f_bus"]
        bus_to = branch["t_bus"]

        if !(data[:bus][bus_fr]["zone"] == data[:bus][bus_to]["zone"])
            push!(coupled, l)
        else
            push!(branches[data[:bus][bus_fr]["zone"]], l)
        end
    end
    
    # Assigning each generator to a region
    gens = Dict{Int, Vector{Int64}}()
    
    for i in zones
        gens[i] = Vector{Int64}()
    end
    for (i, gen) in data[:gen]
        region = data[:bus][gen["gen_bus"]]["zone"]
        push!(gens[region], i)
    end
    
    # Assigning each load to a region
    loads = Dict{Int, Vector{Int64}}()
    
    for i in zones
        loads[i] = Vector{Int64}()
    end
    for (i, load) in data[:load]
        region = data[:bus][load["load_bus"]]["zone"]
        push!(loads[region], i)
    end
    
    # Assigning each shunt to a region
    shunts = Dict{Int, Vector{Int64}}()
    
    for i in zones
        shunts[i] = Vector{Int64}()
    end
    for (i, shunt) in data[:shunt]
        region = data[:bus][shunt["shunt_bus"]]["zone"]
        push!(shunts[region], i)
    end
    
    return (zones, buses, branches, coupled, gens, loads, shunts) 
end

function data_prep(network_data)
    
    (zones, buses, branches, coupled, gens, loads, shunts) = define_regions(data)

    buses_reg = []
    for (k, bus_arr) in buses
        for bus in bus_arr
            push!(buses_reg, (k, bus))
        end
    end
            
    branches_fr  = []
    branches_all = []
    for (k, l_arr) in branches
        for l in l_arr
            i = data[:branch][l]["f_bus"]
            j = data[:branch][l]["t_bus"]
            push!(branches_fr, (k, l, i, j))
            push!(branches_all, (k, l, i, j))        
            push!(branches_all, (k, l, j, i))  
        end
    end

    # branches_fr
    # branches_all

    coupled_fr  = []
    coupled_all = []
    
    extra_buses   = []
    coupled_buses = []

    for l in coupled
        f_bus = data[:branch][l]["f_bus"]
        t_bus = data[:branch][l]["t_bus"]
        f_reg = data[:bus][f_bus]["zone"]
        t_reg = data[:bus][t_bus]["zone"]    

        push!(coupled_fr, [(f_reg, l, f_bus, t_bus), (t_reg, l, f_bus, t_bus)] )
        push!(coupled_all, [(f_reg, l, f_bus, t_bus), (t_reg, l, f_bus, t_bus)] )
        push!(coupled_all, [(f_reg, l, t_bus, f_bus), (t_reg, l, t_bus, f_bus)] )
        
        if !((f_reg, t_bus) in buses_reg)
            push!(extra_buses, (f_reg, t_bus))
            push!(coupled_buses, [(t_reg, t_bus) (f_reg, t_bus)])
        end
        if !((t_reg, f_bus) in buses_reg)
            push!(extra_buses, (t_reg, f_bus))
            push!(coupled_buses, [(t_reg, f_bus) (f_reg, f_bus)])
        end
        
    end
    
#     println(length(extra_buses))
#     println(length(coupled_buses))

    gen_reg = []
    for (k, gen_arr) in gens
        for gen in gen_arr
            push!(gen_reg, (k, gen))
        end
    end

    load_reg = []
    for (k, load_arr) in loads
        for load in load_arr
            push!(load_reg, (k, load))
        end
    end

    shunt_reg = []
    for (k, shunt_arr) in shunts
        for shunt in shunt_arr
            push!(shunt_reg, (k, shunt))
        end
    end
    
    
    
    return Dict("regions" => zones, "branches_fr" => branches_fr, "branches_all" => branches_all, "buses" => buses_reg
                ,"coupled_fr" => coupled_fr, "coupled_all" => coupled_all
                ,"gens" => gen_reg, "loads" => load_reg, "shunts" => shunt_reg
                ,"extra_buses" => extra_buses, "coupled_buses" => coupled_buses)
    
end

function ac_opf_decentralized(data)
    model = Model(Ipopt.Optimizer)
    
    ref   = define_regions(data)
    ref   = data_prep(ref)
    
    # Regions
    regions = ref["regions"]
    
    # Buses
    buses = ref["buses"]
    extra_buses   = ref["extra_buses"]
    coupled_buses = ref["coupled_buses"]
    all_buses     = vcat(buses, extra_buses)
    
    # Intra-Regional Branches
    branches_fr  = ref["branches_fr"]
    branches_all = ref["branches_all"]
    
    # Inter-Regional Branches
    coupled_fr   = ref["coupled_fr"]
    coupled_all  = ref["coupled_all"]
    
    # All Branches
    inter_intra_fr = vcat(collect(Iterators.flatten(coupled_fr)), branches_fr)
    inter_intra    = vcat(collect(Iterators.flatten(coupled_all)), branches_all)
    
    # Other
    gens   = ref["gens"]
    loads  = ref["loads"]
    shunts = ref["shunts"]
    
    # Power Generated
    @variable(model, data[:gen][i]["pmin"] 
        <= pg[(k, i) in gens] 
        <= data[:gen][i]["pmax"])
    
    @variable(model, data[:gen][i]["qmin"] 
        <= qg[(k, i) in gens] 
        <= data[:gen][i]["qmax"])
    
    # Bus Voltage Magnitude and Angles
    @variable(model, data[:bus][i]["vmin"] <= vm[(k, i) in all_buses] <= data[:bus][i]["vmax"])
    @variable(model, va[(k, i) in all_buses])
    
    # Line Power Flow, both intra and inter-regional
    @variable(model, -data[:branch][l]["rate_a"] <= pf[(k, l, i, j) in inter_intra] <= data[:branch][l]["rate_a"])
    @variable(model, -data[:branch][l]["rate_a"] <= qf[(k, l, i, j) in inter_intra] <= data[:branch][l]["rate_a"])
    
    # Objective Function
    @objective(model, Min, 
        sum(data[:gen][i]["cost"][1] * pg[(k, i)]^2 
            + data[:gen][i]["cost"][2] * pg[(k, i)] 
            + data[:gen][i]["cost"][3] for (k, i) in gens)
    )
    
    for ref in keys(data[:ref_buses])
        @constraint(model, va[(data[:ref_buses][ref]["zone"], ref)] == 0.0)
    end
    
    ##### INTRA-REGIONAL CONSTRAINTS #####
    for (k, l, i, j) in branches_fr
        
        branch = data[:branch][l]

        g_fr = branch["g_fr"]
        g_to = branch["g_to"]
        b_fr = branch["b_fr"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2

        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        # Line Power Flow Constraints
        @constraint(model, pf[(k, l, i, j)]^2 + qf[(k, l, i, j)]^2 <= branch["rate_a"]^2)
        @constraint(model, pf[(k, l, j, i)]^2 + qf[(k, l, j, i)]^2 <= branch["rate_a"]^2)
        
        # Angle Differences
        @constraint(model, branch["angmin"] <= va[(k, i)] - va[(k, j)] <= branch["angmax"])

        # Ohm's Law, real, forward
        @NLconstraint(model, pf[(k, l, i, j)] 
                        == (g + g_fr)/tm * vm[(k, i)]^2 
                            + (-g * tr + b * ti)/tm * (vm[(k, i)]*vm[(k, j)]*cos(va[(k, i)]-va[(k, j)])) 
                            + (-b * tr - g * ti)/tm * (vm[(k, i)]*vm[(k, j)]*sin(va[(k, i)]-va[(k, j)]))
        )
        # Ohm's Law real, backward
        @NLconstraint(model, pf[(k, l, j, i)] 
                        == (g + g_to) * vm[(k, j)]^2 
                            + (-g * tr - b * ti)/tm * (vm[(k, i)]*vm[(k, j)]*cos(va[(k, j)]-va[(k, i)])) 
                            + (-b * tr + g * ti)/tm * (vm[(k, i)]*vm[(k, j)]*sin(va[(k, j)]-va[(k, i)]))
        )
        # Ohm's Law Imaginary, forward
        @NLconstraint(model, qf[(k, l, i, j)] 
                        == -(b + b_fr)/tm * vm[(k, i)]^2 
                            - (-b * tr - g * ti)/tm * (vm[(k, i)]*vm[(k, j)]*cos(va[(k, i)]-va[(k, j)])) 
                            + (-g * tr + b * ti)/tm * (vm[(k, i)]*vm[(k, j)]*sin(va[(k, i)]-va[(k, j)]))
        )
        # Ohm's Law Imaginary, backward
        @NLconstraint(model, qf[(k, l, j, i)] 
                        == -(b + b_to) * vm[(k, j)]^2 
                            - (-b * tr + g * ti)/tm * (vm[(k, i)]*vm[(k, j)]*cos(va[(k, j)]-va[(k, i)])) 
                            + (-g * tr - b * ti)/tm * (vm[(k, i)]*vm[(k, j)]*sin(va[(k, j)]-va[(k, i)]))
        )
    end
    
    # Not sure if I can do. Do I need to include duplicate buses?
    for (k, i) in buses
        # Real Power, Kirchoff's Law
        p_gen   = sum(pg[(k, a)] for a in data[:bus_gens][i]; init=0)
        p_load  = sum(data[:load][a]["pd"] for a in data[:bus_loads][i]; init=0)
        p_shunt = sum(data[:shunt][a]["gs"] for a in data[:bus_shunts][i]; init=0) * vm[(k, i)]^2
        p_line  = sum(pf[(k, l, i, j)] for (l, i, j) in data[:bus_arcs][i]; init=0)
        
        @constraint(model, p_gen - p_load - p_shunt == p_line)

        # Imaginary Power, Kirchoff's Law
        q_gen   = sum(qg[(k, a)] for a in data[:bus_gens][i]; init=0)
        q_load  = sum(data[:load][a]["qd"] for a in data[:bus_loads][i]; init=0)
        q_shunt = sum(data[:shunt][a]["bs"] for a in data[:bus_shunts][i]; init=0) * vm[(k, i)]^2
        q_line  = sum(qf[(k, l, i, j)] for (l, i, j) in data[:bus_arcs][i]; init=0)
        
        @constraint(model, q_gen - q_load + q_shunt == q_line)
    end 
    
    #### INTER-REGIONAL CONSTRAINTS ####
    for (k, l, i, j) in collect(Iterators.flatten(coupled_fr))
        
        branch = data[:branch][l]

        g_fr = branch["g_fr"]
        g_to = branch["g_to"]
        b_fr = branch["b_fr"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2

        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        # Line Power Flow Constraints
        @constraint(model, pf[(k, l, i, j)]^2 + qf[(k, l, i, j)]^2 <= branch["rate_a"]^2)
        @constraint(model, pf[(k, l, j, i)]^2 + qf[(k, l, j, i)]^2 <= branch["rate_a"]^2)

        # Angle Differences
        @constraint(model, branch["angmin"] <= va[(k, i)] - va[(k, j)] <= branch["angmax"])

        # Ohm's Law, real, forward
        @NLconstraint(model, pf[(k, l, i, j)] 
                        == (g + g_fr)/tm * vm[(k, i)]^2 
                            + (-g * tr + b * ti)/tm * (vm[(k, i)]*vm[(k, j)]*cos(va[(k, i)]-va[(k, j)])) 
                            + (-b * tr - g * ti)/tm * (vm[(k, i)]*vm[(k, j)]*sin(va[(k, i)]-va[(k, j)]))
        )
        # Ohm's Law real, backward
        @NLconstraint(model, pf[(k, l, j, i)] 
                        == (g + g_to) * vm[(k, j)]^2 
                            + (-g * tr - b * ti)/tm * (vm[(k, i)]*vm[(k, j)]*cos(va[(k, j)]-va[(k, i)])) 
                            + (-b * tr + g * ti)/tm * (vm[(k, i)]*vm[(k, j)]*sin(va[(k, j)]-va[(k, i)]))
        )
        # Ohm's Law Imaginary, forward
        @NLconstraint(model, qf[(k, l, i, j)] 
                        == -(b + b_fr)/tm * vm[(k, i)]^2 
                            - (-b * tr - g * ti)/tm * (vm[(k, i)]*vm[(k, j)]*cos(va[(k, i)]-va[(k, j)])) 
                            + (-g * tr + b * ti)/tm * (vm[(k, i)]*vm[(k, j)]*sin(va[(k, i)]-va[(k, j)]))
        )
        # Ohm's Law Imaginary, backward
        @NLconstraint(model, qf[(k, l, j, i)] 
                        == -(b + b_to) * vm[(k, j)]^2 
                            - (-b * tr + g * ti)/tm * (vm[(k, i)]*vm[(k, j)]*cos(va[(k, j)]-va[(k, i)])) 
                            + (-g * tr - b * ti)/tm * (vm[(k, i)]*vm[(k, j)]*sin(va[(k, j)]-va[(k, i)]))
        )
    end
    
    # Consensus Constraints
    
    # Branches
    for x in coupled_fr
        br_1, br_2 = x
        println(@constraint(model, pf[br_1] == pf[br_2]) )
        println(@constraint(model, qf[br_1] == qf[br_2]) )
    end
    
    # Buses
    for x in coupled_buses
        bus_1, bus_2 = x
        println(@constraint(model, vm[bus_1] == vm[bus_2]) )
        println(@constraint(model, va[bus_1] == va[bus_2]) )
    end
    
    return model
    
end

# Example Use Case
model = ac_opf_decentralized(data)
optimize!(model)
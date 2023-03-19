using PGLib
using PowerModels
using Ipopt
using JuMP
using JSON

# Using the data from PGLib passed into the function
function define_regions(network_data)

    regions = Set(bus_dict["zone"] for (bu, bus_dict) in network_data[:bus])

    # Need to collect inner buses vs outer ones
    inner_bus = Dict()
    outer_bus = Dict()
    frgn_bus  = Dict()
    inner_br  = Dict()
    outer_br  = Dict()
    for k in regions
        inner_bus[k] = Dict()
        outer_bus[k] = Dict()
        frgn_bus[k]  = Dict()
        
        inner_br[k] = Dict()
        outer_br[k] = Dict()
    end

    for (l, br_dict) in network_data[:branch]
        # Track which bus is on each side of the branch
        bus_fr = br_dict["f_bus"]
        bus_to = br_dict["t_bus"]

        # Track which region each of the corresponding buses is in
        reg_fr = network_data[:bus][bus_fr]["zone"]
        reg_to = network_data[:bus][bus_to]["zone"]

        if (reg_fr == reg_to)
            
            # Bus Info
            inner_bus[reg_fr][bus_fr] = network_data[:bus][bus_fr]
            inner_bus[reg_to][bus_to] = network_data[:bus][bus_to]
            
            # Branch Info
            inner_br[reg_fr][l] = network_data[:branch][l]
        else
            # Denote buses which are within the region, on the border of a different one
            outer_bus[reg_fr][bus_fr] = network_data[:bus][bus_fr]
            outer_bus[reg_to][bus_to] = network_data[:bus][bus_to]

            # Denotes buses which are outside of the region, on the border with the current region
            frgn_bus[reg_to][bus_fr] = network_data[:bus][bus_fr]
            frgn_bus[reg_fr][bus_to] = network_data[:bus][bus_to]
            
            # Denotes branches which connect regions
            outer_br[reg_fr][l] = network_data[:branch][l]
            outer_br[reg_to][l] = network_data[:branch][l]
        end

    end
    
    # Add each generator to a region
    gen_dict = Dict()
    for k in regions
        gen_dict[k] = Dict()
    end
    
    for (bus, gens) in network_data[:bus_gens]
        
        region = network_data[:bus][bus]["zone"]
        
        for gen in gens
            gen_dict[region][gen] = network_data[:gen][gen]
        end
        
    end
    
    # Add each load to a region
    load_dict = Dict()
    for k in regions
        load_dict[k] = Dict()
    end
    
    for (bus, loads) in network_data[:bus_loads]
        
        region = network_data[:bus][bus]["zone"]
        
        for load in loads
            load_dict[region][load] = network_data[:load][load]
        end
        
    end

    # Add each shunt to a region
    shunt_dict = Dict()
    for k in regions
        shunt_dict[k] = Dict()
    end
    
    for (bus, shunts) in network_data[:bus_shunts]
        
        region = network_data[:bus][bus]["zone"]
        
        for shunt in shunts
            shunt_dict[region][shunt] = network_data[:shunt][shunt]
        end
        
    end
    
    final_dict = Dict()
    
    for k in regions
        final_dict[k] = Dict("bus"    => Dict("inner" => inner_bus[k], "outer" => outer_bus[k], "foreign" => frgn_bus[k])
                            ,"branch" => Dict("inner" => inner_br[k], "outer" => outer_br[k])
                            ,"gen"    => gen_dict[k]
                            ,"load"   => load_dict[k]
                            ,"shunt"  => shunt_dict[k]                              
                            )
    end

    return final_dict

end

function ac_opf(network_data, ALL_DATA, lmda, consensus, region, rho, prnt_lvl = 0)
    model = Model(Ipopt.Optimizer)

    if prnt_lvl == 0
        set_optimizer_attribute(model, "print_level", prnt_lvl)
    end
    
    regional_data = network_data[region]
    
    # Retrieve all non-inter-regional branches
    inner_br = regional_data["branch"]["inner"]
    
    # Retrieve all inter-regional branches
    outer_br = regional_data["branch"]["outer"]
    
    # Retrieve all buses not connected to an inter-regional branch
    inner_bus = regional_data["bus"]["inner"]
    
    # Retrieve all buses still in the region while also connected to inter-regional branches
    outer_bus = regional_data["bus"]["outer"]
    
    # Find all buses within the region
    local_bus = merge(inner_bus, outer_bus)
    
    # Retrieve all foreign buses
    frgn_bus = regional_data["bus"]["foreign"]
    
    # Retrieve generators, loads, and shunts
    gens   = regional_data["gen"]
    loads  = regional_data["load"]
    shunts = regional_data["shunt"]

    
    # ===
    # Local Variables
    # ===
    
    # Voltage
    @variable(model, vm[keys(local_bus)]) # Magnitude
    @variable(model, va[keys(local_bus)]) # Angle
    
    
#     for ref in keys(z1[:ref_buses])
#         @constraint(model, va[ref] == 0.0)
#     end

    
    # Power Generated
    @variable(model, pg[keys(gens)]) # Active
    @variable(model, qg[keys(gens)]) # Reactive
    
    # Active Line Power Flow
    @variable(model, pf1[keys(inner_br)]) # Forward
    @variable(model, pf2[keys(inner_br)]) # Backward
    
    # Reactive Line Power Flow
    @variable(model, qf1[keys(inner_br)]) # Forward
    @variable(model, qf2[keys(inner_br)]) # Backward
    
    # ===
    # Foreign Variables
    # ===
    
    # Voltage
    @variable(model, vm_f[keys(frgn_bus)]) # Magnitude
    @variable(model, va_f[keys(frgn_bus)]) # Angle
    
    # Active Line Power Flow
    @variable(model, pf1_f[keys(outer_br)]) # Forward
    @variable(model, pf2_f[keys(outer_br)]) # Backward
    
    # Reactive Line Power Flow
    @variable(model, qf1_f[keys(outer_br)]) # Forward
    @variable(model, qf2_f[keys(outer_br)]) # Backward
        
    # Generator Starting Values & Bounds
    for (g, gen_dict) in gens
        JuMP.set_start_value(pg[g], gen_dict["pmin"])
        JuMP.set_start_value(qg[g], gen_dict["qmin"])
        
        # Real Power Bounds
        JuMP.set_lower_bound(pg[g], gen_dict["pmin"])
        JuMP.set_upper_bound(pg[g], gen_dict["pmax"])
        
        # Reactive Power Bounds
        JuMP.set_lower_bound(qg[g], gen_dict["qmin"])
        JuMP.set_upper_bound(qg[g], gen_dict["qmax"])
    end
    
    ###
    # Voltage Starting Value & Bounds
    ###
    
    # Local
    for (bu, bus_dict) in local_bus
        JuMP.set_start_value(vm[bu], 1.0)
        JuMP.set_start_value(va[bu], 0.0)
        
        JuMP.set_lower_bound(vm[bu], bus_dict["vmin"])
        JuMP.set_upper_bound(vm[bu], bus_dict["vmax"])

    end
    # Foreign
    for (bu, bus_dict) in frgn_bus
        JuMP.set_start_value(vm_f[bu], 1.0)
        JuMP.set_start_value(va_f[bu], 0.0)
        
        JuMP.set_lower_bound(vm_f[bu], bus_dict["vmin"])
        JuMP.set_upper_bound(vm_f[bu], bus_dict["vmax"])
    end
    
    # Maximum Angle Differences
    # Local
    for (line, line_dict) in inner_br
        f_bus = line_dict["f_bus"]
        t_bus = line_dict["t_bus"]
        
        @constraint(model, line_dict["angmin"] <= va[f_bus] - va[t_bus] <= line_dict["angmax"])
    end
    
    # Foreign
    for (line, line_dict) in outer_br
        
        f_bus = line_dict["f_bus"]
        t_bus = line_dict["t_bus"]
        
        if (f_bus in keys(frgn_bus))
            @constraint(model, line_dict["angmin"] <= va_f[f_bus] - va[t_bus] <= line_dict["angmax"])
        elseif (t_bus in keys(frgn_bus))
            @constraint(model, line_dict["angmin"] <= va[f_bus] - va_f[t_bus] <= line_dict["angmax"])
        else
            # Sanity Check
            println("buses have not been divided correctly")
        end
    end
    
    # Slack bus
    # May need to add
    
    # Local Thermal Line Limits
    for (line, line_dict) in inner_br
        
        maxCapacity = line_dict["rate_a"]
        
        @NLconstraint(model, pf1[line]^2 + qf1[line]^2 <= maxCapacity^2) # From        
        @NLconstraint(model, pf2[line]^2 + qf2[line]^2 <= maxCapacity^2) # To    
    
    end
    
    # Foreign Thermal Line Limites
    for (line, line_dict) in outer_br
        
        maxCapacity = line_dict["rate_a"]
        
        @NLconstraint(model, pf1_f[line]^2 + qf1_f[line]^2 <= maxCapacity^2) # From        
        @NLconstraint(model, pf2_f[line]^2 + qf2_f[line]^2 <= maxCapacity^2) # To    
        
    end
    
    # ===
    # Ohm's Law
    # ===
    
    # Local Branches
    for (e, line_dict) in inner_br
        g, b   = PowerModels.calc_branch_y(line_dict)
        
        g1 = line_dict["g_fr"]
        g2 = line_dict["g_to"]
        b1 = line_dict["b_fr"]
        b2 = line_dict["b_to"]
        i  = line_dict["f_bus"]
        j  = line_dict["t_bus"]
#         T_m = 1
#         T_R = 1
#         T_I = 0
        T_m = line_dict["tap"]
        T_R, T_I = PowerModels.calc_branch_t(line_dict)

        @NLconstraint(model,
            pf1[e] == (
                (1/T_m^2) * (g + g1) * vm[i]^2
                + ((-g * T_R + b * T_I)/T_m^2) * (vm[i] * vm[j]) * cos(va[i] - va[j])
                + ((-b * T_R - g * T_I)/T_m^2) * (vm[i] * vm[j]) * sin(va[i] - va[j])
            )
        )
        @NLconstraint(model,
            qf1[e] == (
                - (1/T_m^2) * (b + b1) * vm[i]^2
                - ((-b * T_R - g * T_I)/T_m^2) * (vm[i] * vm[j]) * cos(va[i] - va[j])
                + ((-g * T_R + b * T_I)/T_m^2) * (vm[i] * vm[j]) * sin(va[i] - va[j])
            )
        )
        @NLconstraint(model,
            pf2[e] == (
                (g + g2) * vm[j]^2
                + ((-g * T_R - b * T_I)/T_m^2) * (vm[i] * vm[j]) * cos(va[i] - va[j])
                + ((-b * T_R + g * T_I)/T_m^2) * (vm[i] * vm[j]) * sin(-va[i] + va[j])
            )
        )
        @NLconstraint(model,
            qf2[e] == (
                -(b + b2) * vm[j]^2
                - ((-b * T_R + g * T_I)/T_m^2) * (vm[i] * vm[j]) * cos(va[i] - va[j])
                + ((-g * T_R - b * T_I)/T_m^2) * (vm[i] * vm[j]) * sin(-va[i] + va[j])
            )
        )
    end
    
    # Foreign Branches
    for (e, line_dict) in outer_br
        g, b   = PowerModels.calc_branch_y(line_dict)
        
        g1 = line_dict["g_fr"]
        g2 = line_dict["g_to"]
        b1 = line_dict["b_fr"]
        b2 = line_dict["b_to"]
        T_m = line_dict["tap"]
        T_R, T_I = PowerModels.calc_branch_t(line_dict)
        
        i = line_dict["f_bus"]
        j = line_dict["t_bus"]
        
        if (i in keys(frgn_bus))
            vmi = vm_f[i]
            vmj = vm[j]
            vai = va_f[i]
            vaj = va[j]
        elseif (j in keys(frgn_bus))
            vmi = vm[i]
            vmj = vm_f[j]
            vai = va[i]
            vaj = va_f[j]
        else
            println("something is wrong in Ohm's Law!")
        end
        
            
#         T_m = 1
#         T_R = 1
#         T_I = 0


        @NLconstraint(model,
            pf1_f[e] == (
                (1/T_m^2) * (g + g1) * vmi^2
                + ((-g * T_R + b * T_I)/T_m^2) * (vmi * vmj) * cos(vai - vaj)
                + ((-b * T_R - g * T_I)/T_m^2) * (vmi * vmj) * sin(vai - vaj)
            )
        )
        @NLconstraint(model,
            qf1_f[e] == (
                - (1/T_m^2) * (b + b1) * vmi^2
                - ((-b * T_R - g * T_I)/T_m^2) * (vmi * vmj) * cos(vai - vaj)
                + ((-g * T_R + b * T_I)/T_m^2) * (vmi * vmj) * sin(vai - vaj)
            )
        )
        @NLconstraint(model,
            pf2_f[e] == (
                (g + g2) * vmj^2
                + ((-g * T_R - b * T_I)/T_m^2) * (vmi * vmj) * cos(vai - vaj)
                + ((-b * T_R + g * T_I)/T_m^2) * (vmi * vmj) * sin(-vai + vaj)
            )
        )
        @NLconstraint(model,
            qf2_f[e] == (
                -(b + b2) * vmj^2
                - ((-b * T_R + g * T_I)/T_m^2) * (vmi * vmj) * cos(vai - vaj)
                + ((-g * T_R - b * T_I)/T_m^2) * (vmi * vmj) * sin(-vai + vaj)
            )
        )
    end
    
    # =======================
    
    
    # ===
    # Kirchoff's
    # ===
    
    # changed all these to be local buses instead of inner buses
    gen_sum   = Dict(bu => [] for (bu, bu_dict) in local_bus)
    load_sum  = Dict(bu => [] for (bu, bu_dict) in local_bus)
    shunt_sum = Dict(bu => [] for (bu, bu_dict) in local_bus)
    br_in     = Dict(bu => [] for (bu, bu_dict) in local_bus)
    br_out    = Dict(bu => [] for (bu, bu_dict) in local_bus)
    
    br_out_frgn = Dict(bu => [] for (bu, bu_dict) in local_bus)
    br_in_frgn  = Dict(bu => [] for (bu, bu_dict) in local_bus)
    
    # 
    for (g, gen) in gens
        bu = gen["gen_bus"]
        push!(gen_sum[bu], g)
    end
    
    for (l, load) in loads
        bu = load["load_bus"]
        push!(load_sum[bu], l)
    end
    
    for (s, shunt) in shunts
        bu = shunt["shunt_bus"]
        push!(shunt_sum[bu], s)
    end
       
    # 
    for (br, branch) in inner_br
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        
        # Branch leaves from this bus...
        push!(br_out[f_bus], br)
        
        # ...and arrives at this one
        push!(br_in[t_bus], br)
    end
    
    for (br, branch) in outer_br
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        
        # Only need to worry about the local side of the line (ie flow of foreign bus not needed)
        if (f_bus in keys(frgn_bus))
            push!(br_in_frgn[t_bus], br)
        elseif (t_bus in keys(frgn_bus))
            push!(br_out_frgn[f_bus], br)
        # Sanity Check
        else
            println("Issue with Kirchoff's regioning")
        end
        
    end
    
    # This could be a source of issues. Works with one region, but difficult to check with multiple.
    # Loop seems to be working correctly in multi-regional case, but if results don't make sense, this may be a
    #    good place to start
    for (i, bus) in local_bus
        
        # Active Power
        @constraint(model,
            sum(pg[g] for g in gen_sum[i])
            - sum(pf2[e] for e in br_in[i])
            - sum(pf2_f[e] for e in br_in_frgn[i])
            ==
            sum(loads[l]["pd"] for l in load_sum[i])
            + sum(pf1[e] for e in br_out[i])
            + sum(pf1_f[e] for e in br_out_frgn[i])
            + sum(shunts[s]["gs"] for s in shunt_sum[i]) * vm[i]^2
        )

        # Reactive Power
        @constraint(model,
            sum(qg[g] for g in gen_sum[i])
            - sum(qf2[e] for e in br_in[i])
            - sum(qf2_f[e] for e in br_in_frgn[i])
            ==
            sum(loads[l]["qd"] for l in load_sum[i])
            + sum(qf1[e] for e in br_out[i])
            + sum(qf1_f[e] for e in br_out_frgn[i])
            - sum(shunts[s]["bs"] for s in shunt_sum[i]) * vm[i]^2
        )
    end
    
    cst_term = []

    # Cost term values
    for (g, gen) in gens
        push!(cst_term, gen["cost"][1] * pg[g]^2 + gen["cost"][2] * pg[g] + gen["cost"][3])
    end

    
    @objective(model, Min, sum(cst_term; init = 0))


    return model
    
end


function update_obj!(model, network_data, lmda, consensus, region, rho)

    regional_data = network_data[region]
    
    # Retrieve all non-inter-regional branches
    inner_br = regional_data["branch"]["inner"]
    
    # Retrieve all inter-regional branches
    outer_br = regional_data["branch"]["outer"]
    
    # Retrieve all buses not connected to an inter-regional branch
    inner_bus = regional_data["bus"]["inner"]
    
    # Retrieve all buses still in the region while also connected to inter-regional branches
    outer_bus = regional_data["bus"]["outer"]
    
    # Find all buses within the region
    local_bus = merge(inner_bus, outer_bus)
    
    # Retrieve all foreign buses
    frgn_bus = regional_data["bus"]["foreign"]
    
    # Retrieve generators, loads, and shunts
    gens   = regional_data["gen"]
    loads  = regional_data["load"]
    shunts = regional_data["shunt"]

    cst_term = []
    lng_term = []
    pen_term = []
    
    # variables: pf1, pf1_f, pf2, pf2_f, pg, qf1, qf1_f, qf2, qf2_f
    
    # Active power generated
    pg   = model[:pg]
    
    # Bus variables
    vm   = model[:vm]
    vm_f = model[:vm_f]
    va   = model[:va]
    va_f = model[:va_f]
    
    # Branch variables
    pf1_f = model[:pf1_f]
    pf2_f = model[:pf2_f]
    qf1_f = model[:qf1_f]
    qf2_f = model[:qf2_f]
    
    for (g, gen) in gens
        push!(cst_term, gen["cost"][1] * pg[g]^2 + gen["cost"][2] * pg[g] + gen["cost"][3])
    end
    
    # lmda terms
    for (e, line_dict) in outer_br
        f_bus = line_dict["f_bus"]
        t_bus = line_dict["t_bus"]
        
        if (f_bus in keys(frgn_bus))
            vmi = lmda[region]["vm"][f_bus] * vm_f[f_bus]
            vai = lmda[region]["va"][f_bus] * va_f[f_bus]
            
            vmj = lmda[region]["vm"][t_bus] * vm[t_bus]            
            vaj = lmda[region]["va"][t_bus] * va[t_bus]
            
        elseif (t_bus in keys(frgn_bus))
            vmi = lmda[region]["vm"][f_bus] * vm[f_bus]
            vai = lmda[region]["va"][f_bus] * va[f_bus]
            
            vmj = lmda[region]["vm"][t_bus] * vm_f[t_bus]            
            vaj = lmda[region]["va"][t_bus] * va_f[t_bus]
            
        else
            println("lang function is wrong")
        end
        
        push!(lng_term, vmi, vai, vmj, vaj)
        
        pf1 = lmda[region]["pf_1"][e] * pf1_f[e]
        pf2 = lmda[region]["pf_2"][e] * pf2_f[e]
        
        qf1 = lmda[region]["qf_1"][e] * qf1_f[e]
        qf2 = lmda[region]["qf_2"][e] * qf2_f[e]
        
        push!(lng_term, pf1, pf2, qf1, qf2)
        
    end
    
    # consensus terms
    for (e, line_dict) in outer_br
        f_bus = line_dict["f_bus"]
        t_bus = line_dict["t_bus"]
        
        if (f_bus in keys(frgn_bus))
            vmi = (consensus["vm"][f_bus] - vm_f[f_bus])^2
            vai = (consensus["va"][f_bus] - va_f[f_bus])^2
            
            vmj = (consensus["vm"][t_bus] - vm[t_bus])^2
            vaj = (consensus["va"][t_bus] - va[t_bus])^2        
        elseif (t_bus in keys(frgn_bus))
            vmi = (consensus["vm"][f_bus] - vm[f_bus])^2
            vai = (consensus["va"][f_bus] - va[f_bus])^2
            
            vmj = (consensus["vm"][t_bus] - vm_f[t_bus])^2
            vaj = (consensus["va"][t_bus] - va_f[t_bus])^2 
        else
            println("penalty regioning wrong")
        end
        
        push!(pen_term, vmi, vai, vmj, vaj)
        
        pf1 = (consensus["pf_1"][e] - pf1_f[e])^2
        pf2 = (consensus["pf_2"][e] - pf2_f[e])^2
        
        qf1 = (consensus["qf_1"][e] - qf1_f[e])^2
        qf2 = (consensus["qf_1"][e] - qf2_f[e])^2
        
        @objective(model, Min, sum(cst_term; init = 0) + sum(lng_term; init = 0) + rho/2 * sum(pen_term; init = 0))

    end
    
end

function init_consensus(network_data)
    
    regions = [r for (r, data) in network_data]

    consensus = Dict("vm" => Dict(), "va" => Dict()
                        ,"pf_1" => Dict(), "pf_2" => Dict()
                        ,"qf_1" => Dict(), "qf_2" => Dict()
                    )
        
    for k in regions        
        
        lines = network_data[k]["branch"]["outer"]
        
        for (line, line_dict) in lines
            f_bus = line_dict["f_bus"]
            t_bus = line_dict["t_bus"]
            
            # f_bus corresponding consensus
            consensus["va"][f_bus] = 0.0
            consensus["vm"][f_bus] = 1.0

            # t_bus corresponding consensus, okay if duplicates overwrite 
            consensus["va"][t_bus] = 0.0
            consensus["vm"][t_bus] = 1.0
            
            # forward line consensus
            consensus["pf_1"][line] = 0.0
            consensus["qf_1"][line] = 0.0
            
            # backward line consensus
            consensus["pf_2"][line] = 0.0
            consensus["qf_2"][line] = 0.0
            
        end
        
    end
        
    return consensus
    
end

function init_lmda(network_data)
    
    regions = [r for (r, data) in network_data]

    lmda = Dict()
    
    for k in regions
        lmda[k] = Dict("va" => Dict(), "vm" => Dict()
                    ,"pf_1" => Dict(), "pf_2" => Dict()
                    ,"qf_1" => Dict(), "qf_2" => Dict()
                    )
        
        lines = network_data[k]["branch"]["outer"]
        
        for (line, line_dict) in lines
            
            f_bus = line_dict["f_bus"]
            t_bus = line_dict["t_bus"]
            
            # Add lagrangian multiplier for each bus
            lmda[k]["va"][f_bus] = 0.0
            lmda[k]["vm"][f_bus] = 0.0

            # okay if duplicates
            lmda[k]["va"][t_bus] = 0.0
            lmda[k]["vm"][t_bus] = 0.0
            
            # forward line lagrangian multiplier
            lmda[k]["pf_1"][line] = 0.0
            lmda[k]["qf_1"][line] = 0.0
            
            # backward line lagrangian multiplier
            lmda[k]["pf_2"][line] = 0.0
            lmda[k]["qf_2"][line] = 0.0
        end
        
    end
    
    return lmda

end


function update_consensus!(network_data, consensus, model, k)
    
    lines    = network_data[k]["branch"]["outer"]
    frgn_bus = network_data[k]["bus"]["foreign"]

    f_updated = []
    t_updated = []
    
    bus_updated = []
    
    for (line, line_dict) in lines
        f_bus = line_dict["f_bus"]
        t_bus = line_dict["t_bus"]
        
        # We only want to use the consensus value         
        if f_bus in keys(frgn_bus)
            vmi = value(model[:vm_f][f_bus])
            vai = value(model[:va_f][f_bus])

            vmj = value(model[:vm][t_bus])
            vaj = value(model[:va][t_bus])

        elseif t_bus in keys(frgn_bus)
            vmi = value(model[:vm][f_bus])
            vai = value(model[:va][f_bus])

            vmj = value(model[:vm_f][t_bus])
            vaj = value(model[:va_f][t_bus])
        else
            println("issue with districting in consensus")
        end
        
        if !(f_bus in f_updated) & !(f_bus in bus_updated)
        # f_bus corresponding consensus
            consensus["va"][f_bus] = vai
            consensus["vm"][f_bus] = vmi
            
            # Records that bus has been updated. Avoids double updates
            push!(bus_updated, f_bus)
        end

        if !(t_bus in t_updated) & !(t_bus in bus_updated)
        # t_bus corresponding consensus
            consensus["va"][t_bus] = vaj
            consensus["vm"][t_bus] = vmj
            
            # Records that bus has been updated. Avoids double updates
            push!(bus_updated, t_bus)
        end

        # forward line consensus
        consensus["pf_1"][line] = value(model[:pf1_f][line])
        consensus["qf_1"][line] = value(model[:qf1_f][line])

        # backward line consensus
        consensus["pf_2"][line] = value(model[:pf2_f][line])
        consensus["qf_2"][line] = value(model[:qf2_f][line])
        
        # Since multiple lines can connect to a single bus, we need to record
#         push!(already_updated, f_bus, t_bus)    
        push!(f_updated, f_bus)
        push!(t_updated, t_bus)
    end    
end

# This is updating twice for buses which have 2 lines
# Is this the case? Draw this out when I get home and check.
# Need to make sure whichever choice I make is reflected in the objective function later!
function update_lmda!(network_data, lmda, consensus, model, k)
    
    lines    = network_data[k]["branch"]["outer"]
    frgn_bus = network_data[k]["bus"]["foreign"]

    f_updated = []
    t_updated = []
    
    bus_updated = []
    
    for (line, line_dict) in lines
        
        f_bus = line_dict["f_bus"]
        t_bus = line_dict["t_bus"]
        
        # Edge case should be accounted for.
        # if (f_bus in bus_updated) || (t_bus in bus_updated)
        #     println("bus already updated, check edge cases")
        # end
        
        if (f_bus in keys(frgn_bus))
            vmi = value(model[:vm_f][f_bus])
            vai = value(model[:va_f][f_bus])

            vmj = value(model[:vm][t_bus])
            vaj = value(model[:va][t_bus])
            
            push!(bus_updated, f_bus, t_bus)
#             println(vmi)

        elseif (t_bus in keys(frgn_bus))
            vmi = value(model[:vm][f_bus])
            vai = value(model[:va][f_bus])

            vmj = value(model[:vm_f][t_bus])
            vaj = value(model[:va_f][t_bus])
            
            push!(bus_updated, f_bus, t_bus)
#             println(vmi)
            
        # Sanity Check
        else
            println("Districting is wrong in update lmda function")
        end
        
        pf1 = value(model[:pf1_f][line])
        pf2 = value(model[:pf2_f][line])
        qf1 = value(model[:qf1_f][line])
        qf2 = value(model[:qf2_f][line])
        
        if !(f_bus in f_updated)
        # f_bus corresponding lagrangian multiplier
            lmda[k]["va"][f_bus] += (vai - consensus["va"][f_bus])
            lmda[k]["vm"][f_bus] += (vmi - consensus["vm"][f_bus])
        end            

        if !(t_bus in t_updated)
            # t_bus corresponding lagrangian multiplier
            lmda[k]["va"][t_bus] += (vaj - consensus["va"][t_bus])
            lmda[k]["vm"][t_bus] += (vmj - consensus["vm"][t_bus])
        end 
        
        # forward line lagrangian multiplier
        lmda[k]["pf_1"][line] += (pf1 - consensus["pf_1"][line])
        lmda[k]["qf_1"][line] += (qf1 - consensus["qf_1"][line])

        # backward line lagrangian multiplier
        lmda[k]["pf_2"][line] += (pf2 - consensus["pf_2"][line])
        lmda[k]["qf_2"][line] += (qf2 - consensus["qf_2"][line])
        
        # Come back to this
        # if (f_bus in f_updated)
        #     println("multiple lines between 2 regions are connected to a single bus")
        # end
        
        # if (t_bus in t_updated)
        #     println("multiple lines between 2 regions are connected to a single bus")
        # end
        

            
        # Ensures a lambda bus is not updated twice (in case of 2 lines connecting to a region)
        push!(f_updated, f_bus)
        push!(t_updated, t_bus)
        
    end
end

function collect_power_info(network_data, model_dict)
    total_gen  = 0
    total_cost = 0

    regions = keys(network_data)
    
    for k in regions
        
        region_data = network_data[k]
        model = model_dict[k]
        gens  = region_data["gen"]
        
        
        for (g, gen) in gens
            pg = value(model[:pg][g])
            
            total_gen  += pg
            total_cost += gen["cost"][1] * pg^2 + gen["cost"][2] * pg + gen["cost"][3]
        end
        
    end
    
    return total_gen, total_cost
    
end

function collect_solution(network_data, model_dict, consensus)
    
    # Branch, gen, bus
    
    regions  = keys(network_data)
    solution = Dict("branch" => Dict(), "gen" => Dict(), "bus" => Dict())
    
    # branch gets pf, pt, qf, qt
    # bus gets vm, va
    # gen gets pg, qg
    
    for k in regions
        regional_data = network_data[k]
        model = model_dict[k]

        gens     = regional_data["gen"]
        inner_br = regional_data["branch"]["inner"]
        
        for (e, line_dict) in inner_br
            f_bus = line_dict["f_bus"]
            t_bus = line_dict["t_bus"]
            
            pf  = value(model[:pf1][e])
            pt  = value(model[:pf2][e])
            qf  = value(model[:qf1][e])
            qt  = value(model[:qf2][e])
            
            vm1 = value(model[:vm][f_bus])
            vm2 = value(model[:vm][t_bus])
            
            va1 = value(model[:va][f_bus])
            va2 = value(model[:va][t_bus])
                        
            solution["branch"][e]  = Dict("pf" => pf, "pt" => pt, "qf" => qf, "qt" => qt)
            solution["bus"][f_bus] = Dict("va" => va1, "vm" => vm1)
            solution["bus"][t_bus] = Dict("va" => va2, "vm" => va2)
        end
        
        for (g, gen) in gens
            pg = value(model[:pg][g])
            qg = value(model[:qg][g])
            solution["gen"][g] = Dict("pg" => pg, "qg" => qg)
        end
    end
    
    va_con  = consensus["va"]
    vm_con  = consensus["vm"]
    pf1_con = consensus["pf_1"]
    pf2_con = consensus["pf_2"]
    qf1_con = consensus["qf_1"]
    qf2_con = consensus["qf_2"]
    
#     bus_info = Dict()
#     br_info  = Dict()
    
    for (bus_num, va) in va_con
        solution["bus"][bus_num]["va"] = va
    end
    
    for (bus_num, vm) in vm_con
        solution["bus"][bus_num]["vm"] = vm
    end
    
    for (br_num, pf1) in pf1_con
        if !(br_num in keys(solution["branch"]))
            solution["branch"][br_num] = Dict()
        end
        solution["branch"][br_num]["pf"] = pf1
    end
    
    for (br_num, pf2) in pf2_con
        solution["branch"][br_num]["pt"] = pf2
    end
    
    for (br_num, qf1) in qf1_con
        solution["branch"][br_num]["qf"] = qf1
    end
    
    for (br_num, qf2) in qf2_con
        solution["branch"][br_num]["qt"] = qf2
    end

    return solution
    
end

function retrieve_inf_primal(network_data, model_dict, consensus)
    
    regions  = keys(network_data)
    inter_br = Dict()
    
    inf = []
    
    # Collect all branches
    for k in regions
        
        region_data = network_data[k]
        frgn_bus    = region_data["bus"]["foreign"]
        
        for (e, line_dict) in region_data["branch"]["outer"]
            f_bus = line_dict["f_bus"]
            t_bus = line_dict["t_bus"]
            
            if (f_bus in keys(frgn_bus))
                reg1 = frgn_bus[f_bus]["zone"]
                reg2 = k
                
                # Active Power
                pf1 = abs(value(model_dict[reg1][:pf1_f][e]) - value(model_dict[reg2][:pf1_f][e]))
                pf2 = abs(value(model_dict[reg1][:pf2_f][e]) - value(model_dict[reg2][:pf2_f][e]))
                
                # Reactive Power
                qf1 = abs(value(model_dict[reg1][:qf1_f][e]) - value(model_dict[reg2][:qf1_f][e]))
                qf2 = abs(value(model_dict[reg1][:qf2_f][e]) - value(model_dict[reg2][:qf2_f][e]))
                
                # Starting bus (f_bus)
                vm1 = abs(value(model_dict[reg1][:vm][f_bus]) - value(model_dict[reg2][:vm_f][f_bus]))
                va1 = abs(value(model_dict[reg1][:va][f_bus]) - value(model_dict[reg2][:va_f][f_bus]))
                
                # Ending bus (t_bus)
                vm2 = abs(value(model_dict[reg1][:vm_f][t_bus]) - value(model_dict[reg2][:vm][t_bus]))
                va2 = abs(value(model_dict[reg1][:va_f][t_bus]) - value(model_dict[reg2][:va][t_bus]))
                
                
            elseif (t_bus in keys(frgn_bus))
                reg1 = k
                reg2 = frgn_bus[t_bus]["zone"]
                
                # Active Power
                pf1 = abs(value(model_dict[reg1][:pf1_f][e]) - value(model_dict[reg2][:pf1_f][e]))
                pf2 = abs(value(model_dict[reg1][:pf2_f][e]) - value(model_dict[reg2][:pf2_f][e]))
                
                # Reactive Power
                qf1 = abs(value(model_dict[reg1][:qf1_f][e]) - value(model_dict[reg2][:qf1_f][e]))
                qf2 = abs(value(model_dict[reg1][:qf2_f][e]) - value(model_dict[reg2][:qf2_f][e]))
                
                # Starting bus (f_bus)
                vm1 = abs(value(model_dict[reg1][:vm][f_bus]) - value(model_dict[reg2][:vm_f][f_bus]))
                va1 = abs(value(model_dict[reg1][:va][f_bus]) - value(model_dict[reg2][:va_f][f_bus]))
                
                # Ending bus (t_bus)
                vm2 = abs(value(model_dict[reg1][:vm_f][t_bus]) - value(model_dict[reg2][:vm][t_bus]))
                va2 = abs(value(model_dict[reg1][:va_f][t_bus]) - value(model_dict[reg2][:va][t_bus]))
                
            else
                println("inf primal districting incorrect")
            end
            
            push!(inf, max(pf1, pf2, qf1, qf2, vm1, vm2, va1, va2))
            
        end
        
    end
    
    return maximum(inf)
    
end

function retrieve_inf_dual(current_sol, prev_sol, rho)
    difs = []
    
#     Dict("branch" => Dict(), "gen" => Dict(), "bus" => Dict())
    for (br, new_dict) in current_sol["branch"]
        old_dict = prev_sol["branch"][br]
        
        pf = abs(new_dict["pf"] - old_dict["pf"])
        pt = abs(new_dict["pt"] - old_dict["pt"])
        qf = abs(new_dict["qf"] - old_dict["qf"])
        qt = abs(new_dict["qt"] - old_dict["qt"])
        
        push!(difs, pf, pt, qf, qt)
        
    end
    
    for (bu, new_dict) in current_sol["bus"]
        old_dict = prev_sol["bus"][bu]
        
        va = abs(new_dict["va"] - old_dict["va"])
        vm = abs(new_dict["vm"] - old_dict["vm"])

        push!(difs, va, vm)
    end
    
    for (gen, new_dict) in current_sol["gen"]
        old_dict = prev_sol["gen"][gen]
        
        pg = abs(new_dict["pg"] - old_dict["pg"])
        qg = abs(new_dict["qg"] - old_dict["qg"])

        push!(difs, pg, qg)
    end
    
    return maximum(difs) * rho
    
end

function main_routine(raw_data, rho, maxIter)

#     define_regions(network_data)
    network_data  = define_regions(raw_data)
    regions       = keys(network_data)
    final_results = Dict()
#     init_consensus(network_data)
    consensus = init_consensus(network_data)
    
#     init_lmda(network_data)
    lmda = init_lmda(network_data)
    
    
    for iter in range(start = 1, stop = maxIter)
        
        println(iter)
        iter_results     = Dict()
        final_results[iter] = Dict()
        
        for k in regions
    #         ac_opf(network_data, ALL_DATA, lmda, consensus, region, rho, prnt_lvl = 0)
            model = ac_opf(network_data, raw_data, lmda, consensus, k, rho)
            
#             update_obj!(model, network_data, lmda, consensus, region, rho)
            update_obj!(model, network_data, lmda, consensus, k, rho)
            
            optimize!(model)

            #     update_lmda!(network_data, lmda, consensus, model, k)
            update_lmda!(network_data, lmda, consensus, model, k)

            #     update_consensus!(network_data, consensus, model, k)
            update_consensus!(network_data, consensus, model, k)

            
            iter_results[k] = model
            
        end
        
        # Record current model results
        total_gen, total_cost = collect_power_info(network_data, iter_results)
        primal_inf = retrieve_inf_primal(network_data, iter_results, consensus)
        curr_sol   = collect_solution(network_data, iter_results, consensus)
        
        final_results[iter]["gen"]      = total_gen
        final_results[iter]["cost"]     = total_cost
        final_results[iter]["primal"]   = primal_inf
        final_results[iter]["solution"] = curr_sol
        
        if iter == 1
            final_results[iter]["dual"] = 0
        else
            prev_sol = final_results[iter - 1]["solution"]
            dual_inf = retrieve_inf_dual(curr_sol, prev_sol, rho)
            final_results[iter]["dual"] = dual_inf
        end    
        

    end
    
    return final_results
end

# Takes in a PGLib dictionary and returns the final results
# Parameters: PGLib Dict, rho, max_iterations
function solve_admm_ac_opf(pg_lib_dict, rho, maxIter, export_res = false, export_loc = "")

    PowerModels.standardize_cost_terms!(pg_lib_dict, order=2)
    
    # Adds reasonable rate_a values to branches without them
    PowerModels.calc_thermal_limits!(pg_lib_dict)

    # use build_ref to filter out inactive components
    network_data = PowerModels.build_ref(pg_lib_dict)[:it][:pm][:nw][0]

    results = main_routine(network_data, rho, maxIter)

    if export_res == true
        open(export_loc,"w") do f
            JSON.print(f, results, 4)
        end
    end

    final_result = results[maxIter]

    return Dict("optimizer" => "Ipopt", "objective" => final_result["cost"]
                ,"solution" => final_result["solution"]
                ,"Primal Inf" => final_result["primal"]
                ,"Dual Inf" => final_result["dual"])

end

# Example

# data    = pglib("pglib_opf_case793_goc.m")
# rho1    = 5000
# rho2    = 1
# maxIter = 500

# res = solve_admm_ac_opf(data, rho1, maxIter, true, "case_793_5000.json")
# res = solve_admm_ac_opf(data, rho2, maxIter, true, "case_793_0001.json")
# res = solve_admm_ac_opf(data, 1 * 10^9, 500, true, "case_793_BIG_500.json")
# res

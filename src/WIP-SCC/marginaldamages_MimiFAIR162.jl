
# A component for an emissions pulse to be used in social cost calculations. Computes the `output` vector by adding
#   `add` to `input`. This is similar to the Mimi.adder component, except that it allows missing values to be passed through.

@defcomp emissionspulse begin
    add    = Parameter(index=[time])
    input  = Parameter(index=[time])
    output = Variable(index=[time])

    function run_timestep(p, v, d, t)
        v.output[t] = Mimi.@allow_missing(p.input[t]) + p.add[t]
    end
end

function create_marginal_fundfair16_model(m; year::Int64, pulse_size = 1e7, gas::Symbol = :CO2)

    if gas == :CO2
        
        mm = Mimi.create_marginal_model(m, pulse_size)

        # Add additional emissions to mm.modified
        add_comp!(mm.modified, emissionspulse, after = :emissions)
        nyears = length(Mimi.time_labels(m))
        addem = zeros(nyears) 
        baseyear = 1750
        # if year != nothing 
        #     # pulse is spread over ten years, and emissions components is in Mt so divide by 1e7, and convert from CO2 to C if gas==:CO2 because emissions component is in MtC
            addem[(year - baseyear + 1):(year - baseyear + 1) + 9] .= pulse_size / 1e7 * (gas == :CO2 ? 12/44 : 1)
            # addem[(year - baseyear + 1):(year - baseyear + 1)] .= pulse_size / 1e6 * (gas == :CO2 ? 12/44 : 1) # 1-year pulse
        # end
        set_param!(mm.modified, :emissionspulse, :add, addem)

        # Reconnect the appropriate emissions in m
        connect_param!(mm.modified, :emissionspulse, :input, :emissions, :mco2, zeros(nyears), backup_offset = 1)
        connect_param!(mm.modified, :multiplier, :input, :emissionspulse, :output)
    elseif gas == :CH4
        
        mm = Mimi.create_marginal_model(m, pulse_size/10) # only a 1 year pulse

        pulse_year_index = findall((in)([year]), collect(1750:2300))

        # perturb CH4 emissions
        new_emissions = deepcopy(mm.base[:ch4_cycle, :fossil_emiss_CH₄])
        new_emissions[pulse_year_index] = new_emissions[pulse_year_index] .+ (pulse_size/1e7) # 1MtCH4 pulse

        # update emissions parameter
        Mimi.set_param!(mm.modified, :fossil_emiss_CH₄, new_emissions)

    elseif gas == :N2O

        mm = Mimi.create_marginal_model(m, pulse_size/10) # 1 year pulse

        pulse_year_index = findall((in)([year]), collect(1750:2300))

        # perturb N2O emissions
        new_emissions = deepcopy(mm.base[:n2o_cycle, :fossil_emiss_N₂O])
        new_emissions[pulse_year_index] = new_emissions[pulse_year_index] .+ (pulse_size/1e7 * 28/44) # N2O emissions are in MtN2, multiply by 28/44 to make it a 1MtN2O pulse

        # update emissions parameter
        Mimi.set_param!(mm.modified, :fossil_emiss_N₂O, new_emissions)
        
    end
    
    run(mm)

    return(mm)

end


function compute_fundfair16_sc(mm; rate::Float64, year::Int64, start_year::Int64 = 1750, end_year::Int64 = 2300)
    
    marginaldamage = mm[:impactaggregation, :loss]

    # calculate discount factors
    new_years = collect(start_year:1:end_year)
    year_index = findfirst(isequal(year), new_years)

    df = zeros(length(new_years), 16)
    for i in 1:length(new_years)
        if i >= year_index
            df[i,:] .= 1/(1+rate)^(i-year_index)
        end
    end

    # calculate SCGHG
    sc = sum(skipmissing(marginaldamage .* df))
    return(sc)

end
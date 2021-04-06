using Mimi
using DataFrames

function pad_parameters(data::Array, time_len::Number, begin_padding::Number, end_padding::Number)
    size(data, 1) != time_len ? error("time dimension must match rows") : nothing

    new_data = deepcopy(data)
    if begin_padding != 0
        begin_padding_rows = Array{Union{Missing, Number}}(missing, begin_padding, size(data, 2)) 
        new_data = vcat(begin_padding_rows, new_data)
    end

    if end_padding != 0
        end_padding_rows = Array{Union{Missing, Number}}(missing, end_padding, size(data, 2)) 
        new_data = vcat(new_data, end_padding_rows)
    end

    ndims(data) == 1 ? new_data = new_data[:,1] : nothing 

    return new_data
end

function pad_parameters(data::DataFrame, time_len::Number, begin_padding::Number, end_padding::Number)
    fields = names(data)
    arr = convert(Array, data)
    new_data = pad_parameters(arr, time_len, begin_padding, end_padding) |> DataFrame
    rename!(new_data, Symbol.(fields))
end

# this approach works if we don't want to pull in all of the FUND data, 
# but it's ugly ... should think this through
function update_MimiFUND_params!(m)

    model_first = m.md.first
    model_last = m.md.last

    comp_first = 1950
    comp_last = 2500
    
    data_first = 1950
    data_last = 3000

    # the original FUND parameters end in 3000 so we need to trim off the last 
    # 500 years 
    begin_trim = comp_first - data_first;
    end_trim = data_last - comp_last;

    # the original FUND starts in 1950 and ends in 2500 so we need to pad the 
    # beginning to account for the earlier model start
    begin_padding = comp_first - model_first
    end_padding = model_last - comp_last
    
    for (name, value) in m.md.external_params
        if ~(value isa Mimi.ScalarModelParameter) && (:time in value.dim_names) 
            
            # get data
            data = value.values.data

            # first trim
            if ndims(data) == 1
                data = data[begin_trim + 1:end-end_trim]
            elseif ndims(data) == 2
                data = data[begin_trim + 1:end-end_trim,:]
            end

            # now pad
            padded_data = pad_parameters(data, length(comp_first:comp_last), begin_padding, end_padding)

            # now update parameter
            update_param!(m, name, padded_data)
        end
    end
end

function set_MimiFAIR_params!(m; rcp_scenario::String="RCP85", first::Int=1765, last::Int=2500, F2x::Float64=3.71, TCR::Float64=1.6, ECS::Float64=2.75, d::Array{Float64,1}=[239.0, 4.1])
    
    model_first = m.md.first
    model_last = m.md.last
    
    # ---------------------------------------------
    # Set up some data.
    # ---------------------------------------------

    # Load RCP and other data needed to construct FAIR.
    rcp_emissions, volcano_forcing, solar_forcing, gas_data, gas_fractions, conversions = MimiFAIR.load_fair_data(first, last, rcp_scenario)

    # Calculate CO₂ radiative forcing scaling term (to keep forcing from 2x CO₂ consistent).
    scale_co2_forcing = MimiFAIR.co2_rf_scale(F2x, gas_data[gas_data.gas .== "CO2", :pi_conc][1], gas_data[gas_data.gas .== "N2O", :pi_conc][1])

    # Calcualte coefficients of slow and fast temperature change in each timestep based on user=specified parameters.
    q = MimiFAIR.calculate_q(TCR, ECS, d, F2x)

    # Number of time periods to run model.
    comp_time_len = length(first:last)

    # ---------------------------------------------
    # Padding Data with Time Dimension
    # ---------------------------------------------

    begin_padding = first - model_first
    end_padding = model_last - last

    # pad parameters that have a time dimension
    rcp_emissions = pad_parameters(rcp_emissions, comp_time_len, begin_padding, end_padding)
    volcano_forcing = pad_parameters(volcano_forcing, comp_time_len, begin_padding, end_padding)
    solar_forcing= pad_parameters(solar_forcing, comp_time_len, begin_padding, end_padding)
    gas_fractions = pad_parameters(gas_fractions, comp_time_len, begin_padding, end_padding)

    # ---------------------------------------------
    # Set component parameters
    # ---------------------------------------------

    # ---- Methane Cycle ---- #
    set_param!(m, :MimiFAIR_composite, :fossil_emiss_CH₄, rcp_emissions.CH4)
    set_param!(m, :MimiFAIR_composite, :natural_emiss_CH₄, rcp_emissions.NaturalCH4)
    set_param!(m, :MimiFAIR_composite, :τ_CH₄, 9.3)
    set_param!(m, :MimiFAIR_composite, :fossil_frac, gas_fractions.ch4_fossil)
    set_param!(m, :MimiFAIR_composite, :oxidation_frac, 0.61)
    set_param!(m, :MimiFAIR_composite, :mol_weight_CH₄, gas_data[gas_data.gas .== "CH4", :mol_weight][1])
    set_param!(m, :MimiFAIR_composite, :mol_weight_C, gas_data[gas_data.gas .== "C", :mol_weight][1])
    set_param!(m, :MimiFAIR_composite, :emiss2conc_ch4, conversions[conversions.gases .== "CH4", :emiss2conc][1])

    # ---- Nitrous Oxide Cycle ---- #
    set_param!(m, :MimiFAIR_composite, :fossil_emiss_N₂O, rcp_emissions.N2O)
    set_param!(m, :MimiFAIR_composite, :natural_emiss_N₂O, rcp_emissions.NaturalN2O)
    set_param!(m, :MimiFAIR_composite, :τ_N₂O, 121.0)
    set_param!(m, :MimiFAIR_composite, :emiss2conc_n2o, conversions[conversions.gases .== "N2O", :emiss2conc][1])

    # ---- Carbon Cycle ---- #
    set_param!(m, :MimiFAIR_composite,  :CO2_0, gas_data[gas_data.gas .== "CO2", :pi_conc][1])
    set_param!(m, :MimiFAIR_composite,  :r0, 35.0)
    set_param!(m, :MimiFAIR_composite,  :rC, 0.019)
    set_param!(m, :MimiFAIR_composite,  :rT, 4.165)
    set_param!(m, :MimiFAIR_composite,  :iIRF_max, 97.0)
    set_param!(m, :MimiFAIR_composite,  :a, [0.2173, 0.2240, 0.2824, 0.2763])
    set_param!(m, :MimiFAIR_composite,  :τ_CO₂, [10.0^6, 394.4, 36.54, 4.304])
    set_param!(m, :MimiFAIR_composite,  :E_CO₂, rcp_emissions.FossilCO2 .+ rcp_emissions.OtherCO2)

    # ---- Other Well-Mixed Greenhouse Gas Cycles ---- #
    set_param!(m, :MimiFAIR_composite,  :τ_other_ghg, gas_data[findall((in)(other_ghg_names), gas_data.gas), :lifetimes])    
    set_param!(m, :MimiFAIR_composite,  :emiss_other_ghg, convert(Matrix, rcp_emissions[!,Symbol.(other_ghg_names)]))
    set_param!(m, :MimiFAIR_composite,  :emiss2conc_other_ghg, conversions[findall((in)(other_ghg_names), conversions.gases), :emiss2conc])

    # ---- Methane Radiative Forcing ---- #
    set_param!(m, :MimiFAIR_composite,  :H₂O_share, 0.12)
    set_param!(m, :MimiFAIR_composite,  :scale_CH₄, 1.0)
    set_param!(m, :MimiFAIR_composite,  :a₃, -1.3e-6)
    set_param!(m, :MimiFAIR_composite,  :b₃, -8.2e-6)

    # ---- Nitrous Oxide Radiative Forcing ---- #
    set_param!(m, :MimiFAIR_composite,  :a₂, -8.0e-6)
    set_param!(m, :MimiFAIR_composite,  :b₂, 4.2e-6)
    set_param!(m, :MimiFAIR_composite,  :c₂, -4.9e-6)    

    # ---- Carbon Dioxide Radiative Forcing ---- #
    set_param!(m, :MimiFAIR_composite,  :a₁, -2.4e-7)
    set_param!(m, :MimiFAIR_composite,  :b₁, 7.2e-4)
    set_param!(m, :MimiFAIR_composite,  :c₁, -2.1e-4)
    set_param!(m, :MimiFAIR_composite,  :rf_scale_CO₂, scale_co2_forcing)

    # ---- Other Well-Mixed Greenhouse Gas Radiative Forcings ---- #
    set_param!(m, :MimiFAIR_composite,  :radiative_efficiency, gas_data[findall((in)(other_ghg_names), gas_data.gas), :rad_eff])

    # ---- Tropospheric Ozone Radiative Forcing ---- #
    set_param!(m, :MimiFAIR_composite,  :NOx_emissions, rcp_emissions.NOx)
    set_param!(m, :MimiFAIR_composite,  :CO_emissions, rcp_emissions.CO)
    set_param!(m, :MimiFAIR_composite,  :NMVOC_emissions, rcp_emissions.NMVOC)    
    set_param!(m, :MimiFAIR_composite,  :mol_weight_NO, gas_data[gas_data.gas .== "NO", :mol_weight][1])
    set_param!(m, :MimiFAIR_composite,  :T0, 0.0)    

    # ---- Stratospheric Ozone Radiative Forcing ---- #
    set_param!(m, :MimiFAIR_composite,  :Br, gas_data[findall((in)(ods_names), gas_data.gas), :br_atoms])
    set_param!(m, :MimiFAIR_composite,  :Cl, gas_data[findall((in)(ods_names), gas_data.gas), :cl_atoms])
    set_param!(m, :MimiFAIR_composite,  :FC, gas_data[findall((in)(ods_names), gas_data.gas), :strat_frac])
    set_param!(m, :MimiFAIR_composite,  :δ1, -1.46030698e-5)
    set_param!(m, :MimiFAIR_composite,  :δ2, 2.05401270e-3)
    set_param!(m, :MimiFAIR_composite,  :δ3, 1.03143308)
    set_param!(m, :MimiFAIR_composite,  :ODS₀, gas_data[findall((in)(ods_names), gas_data.gas), :pi_conc])

    # ---- Aerosol Direct Radiative Forcing ---- #
    set_param!(m, :MimiFAIR_composite,  :β_SOx, -6.2227e-3)
    set_param!(m, :MimiFAIR_composite,  :β_CO, 0.0)
    set_param!(m, :MimiFAIR_composite,  :β_NMVOC, -3.8392e-4)
    set_param!(m, :MimiFAIR_composite,  :β_NOx, -1.16551e-3)
    set_param!(m, :MimiFAIR_composite,  :β_BC, 1.601537e-2)
    set_param!(m, :MimiFAIR_composite,  :β_OC, -1.45339e-3)
    set_param!(m, :MimiFAIR_composite,  :β_NH3, -1.55605e-3)
    set_param!(m, :MimiFAIR_composite,  :rf_scale_aero_direct, 0.0)    
    set_param!(m, :MimiFAIR_composite,  :CO_emiss, rcp_emissions.CO)
    set_param!(m, :MimiFAIR_composite,  :NMVOC_emiss, rcp_emissions.NMVOC)
    set_param!(m, :MimiFAIR_composite,  :NH3_emiss, rcp_emissions.NH3)

    # ---- Aerosol Indirect Radiative Forcing ---- #
    set_param!(m, :MimiFAIR_composite,  :ϕ, -1.95011431)
    set_param!(m, :MimiFAIR_composite,  :b_SOx, 0.01107147)
    set_param!(m, :MimiFAIR_composite,  :b_POM, 0.01387492)
    set_param!(m, :MimiFAIR_composite,  :rf_scale_aero_indirect, 0.0)
    model_years = pad_parameters(collect(first:last), comp_time_len, begin_padding, end_padding)
    set_param!(m, :MimiFAIR_composite,  :model_years, model_years)
    set_param!(m, :MimiFAIR_composite,  :SOx_emiss_1765, 1.0)
    set_param!(m, :MimiFAIR_composite,  :BC_OC_emiss_1765, 11.2)
    set_param!(m, :MimiFAIR_composite,  :scale_AR5, true)
    set_param!(m, :MimiFAIR_composite,  :F_1765, -0.3002836449793625)
    set_param!(m, :MimiFAIR_composite,  :F_2011, -1.5236182344467388)

    # ---- Black Carbon on Snow Radiative Forcing ---- #

    # ---- Land Use Change Radiative Forcing ---- #
    set_param!(m, :MimiFAIR_composite,  :landuse_emiss, rcp_emissions.OtherCO2)

    # ---- Contrails Radiative Forcing ---- #    
    set_param!(m, :MimiFAIR_composite,  :frac, gas_fractions.nox_aviation)
    set_param!(m, :MimiFAIR_composite,  :E_ref, 2.946)
    set_param!(m, :MimiFAIR_composite,  :F_ref, 0.0448)
    set_param!(m, :MimiFAIR_composite,  :ref_is_NO2, true)
    set_param!(m, :MimiFAIR_composite,  :mol_weight_NO₂, gas_data[gas_data.gas .== "NO2", :mol_weight][1])    

    # ---- Total Radiative Forcing ---- #
    set_param!(m, :MimiFAIR_composite, :F_volcanic, volcano_forcing)
    set_param!(m, :MimiFAIR_composite, :F_solar, solar_forcing)
    set_param!(m, :MimiFAIR_composite, :F_exogenous, zeros(m.md.last - m.md.first + 1))
    set_param!(m, :MimiFAIR_composite, :efficacy_CO₂, 1.0)
    set_param!(m, :MimiFAIR_composite, :efficacy_CH₄, 1.0)
    set_param!(m, :MimiFAIR_composite, :efficacy_CH₄_H₂O, 1.0)
    set_param!(m, :MimiFAIR_composite, :efficacy_N₂O, 1.0)
    set_param!(m, :MimiFAIR_composite, :efficacy_other_ghg, ones(length(other_ghg_names)))
    set_param!(m, :MimiFAIR_composite, :efficacy_trop_O₃, 1.0)
    set_param!(m, :MimiFAIR_composite, :efficacy_strat_O₃, 1.0)
    set_param!(m, :MimiFAIR_composite, :efficacy_aerosol_direct, 1.0)
    set_param!(m, :MimiFAIR_composite, :efficacy_aerosol_indirect, 1.0)
    set_param!(m, :MimiFAIR_composite, :efficacy_bcsnow, 3.0)
    set_param!(m, :MimiFAIR_composite, :efficacy_landuse, 1.0)
    set_param!(m, :MimiFAIR_composite, :efficacy_contrails, 0.0) # Note: Efficacy set to 0.0 to match default settings in Python version of FAIR.

    # ---- Global Temperature Anomaly ---- #
    set_param!(m, :MimiFAIR_composite,  :d, d)
    set_param!(m, :MimiFAIR_composite,  :q, q)
    set_param!(m, :MimiFAIR_composite,  :F2x, F2x)

    # --- Set common parameters ---- #
    set_param!(m, :CO₂_0, gas_data[gas_data.gas .== "CO2", :pi_conc][1])
    set_param!(m, :CH₄_0, gas_data[gas_data.gas .== "CH4", :pi_conc][1])
    set_param!(m, :N₂O_0, gas_data[gas_data.gas .== "N2O", :pi_conc][1])
    set_param!(m, :other_ghg_0, gas_data[findall((in)(other_ghg_names), gas_data.gas), :pi_conc])
    set_param!(m, :SOx_emiss, rcp_emissions.SOx)
    set_param!(m, :BC_emiss, rcp_emissions.BC)
    set_param!(m, :OC_emiss, rcp_emissions.OC)
    set_param!(m, :NOx_emiss, rcp_emissions.NOx)
    set_param!(m, :fix_pre1850_RCP, true)
    set_param!(m, :mol_weight_N, gas_data[gas_data.gas .== "N", :mol_weight][1])
    set_param!(m, :gtc2ppm, conversions[conversions.gases .== "CO2", :emiss2conc][1])

    return m
end

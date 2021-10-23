using Mimi
using DataFrames
using CSVFiles

"""
    update_MimiFAIR_params!(m; 
                            rcp_scenario::String="RCP85", 
                            start_year::Int=1765, 
                            end_year::Int=2500, 
                            F2x::Float64=3.71, 
                            TCR::Float64=1.6, 
                            ECS::Float64=2.75, 
                            d::Array{Float64,1}=[239.0, 4.1])

Update all parameter settings in the model `m` with default FAIR settings as pulled
from the MimiFAIR repository.  Optional arguments include `rcp_scenario` to indicate
what scenario to use, `start_year` and `end_year` to indicate the model time dimension,
and three parameter settings for the paramters `F2x`, `TCR`, and `ECS`.  See the
MimiFAIR model documentatoin for more.
"""
function update_MimiFAIR_params!(m; rcp_scenario::String="RCP85", start_year::Int=1765, end_year::Int=2500, F2x::Float64=3.71, TCR::Float64=1.6, ECS::Float64=2.75, d::Array{Float64,1}=[239.0, 4.1])
    
    # ---------------------------------------------
    # Set up some data.
    # ---------------------------------------------

    # Load RCP and other data needed to construct FAIR.
    rcp_emissions, volcano_forcing, solar_forcing, gas_data, gas_fractions, conversions = MimiFAIR.load_fair_data(start_year, end_year, rcp_scenario)

    # Names of minor greenhouse gases and ozone-depleting substances.
    other_ghg_names = ["CF4", "C2F6", "C6F14", "HFC23", "HFC32", "HFC43_10", "HFC125", "HFC134a", "HFC143a", "HFC227ea", "HFC245fa", "SF6", "CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]
    ods_names = ["CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]

    # Calculate CO₂ radiative forcing scaling term (to keep forcing from 2x CO₂ consistent).
    scale_co2_forcing = MimiFAIR.co2_rf_scale(F2x, gas_data[gas_data.gas .== "CO2", :pi_conc][1], gas_data[gas_data.gas .== "N2O", :pi_conc][1])

    # Calcualte coefficients of slow and fast temperature change in each timestep based on user=specified parameters.
    q = MimiFAIR.calculate_q(TCR, ECS, d, F2x)

    # Number of time periods to run model.
    n_steps = length(start_year:end_year)

    # ---------------------------------------------
    # Set component parameters
    # ---------------------------------------------

    # ---- Methane Cycle ---- #
    update_param!(m, :ch4_cycle, :fossil_emiss_CH₄, rcp_emissions.CH4)
    update_param!(m, :ch4_cycle, :natural_emiss_CH₄, rcp_emissions.NaturalCH4)
    update_param!(m, :ch4_cycle, :τ_CH₄, 9.3)
    update_param!(m, :ch4_cycle, :fossil_frac, gas_fractions.ch4_fossil)
    update_param!(m, :ch4_cycle, :oxidation_frac, 0.61)
    update_param!(m, :ch4_cycle, :mol_weight_CH₄, gas_data[gas_data.gas .== "CH4", :mol_weight][1])
    update_param!(m, :ch4_cycle, :mol_weight_C, gas_data[gas_data.gas .== "C", :mol_weight][1])
    update_param!(m, :ch4_cycle, :emiss2conc_ch4, conversions[conversions.gases .== "CH4", :emiss2conc][1])

    # ---- Nitrous Oxide Cycle ---- #
    update_param!(m, :n2o_cycle, :fossil_emiss_N₂O, rcp_emissions.N2O)
    update_param!(m, :n2o_cycle, :natural_emiss_N₂O, rcp_emissions.NaturalN2O)
    update_param!(m, :n2o_cycle, :τ_N₂O, 121.0)
    update_param!(m, :n2o_cycle, :emiss2conc_n2o, conversions[conversions.gases .== "N2O", :emiss2conc][1])

    # ---- Carbon Cycle ---- #
    update_param!(m, :co2_cycle, :CO2_0, gas_data[gas_data.gas .== "CO2", :pi_conc][1])
    update_param!(m, :co2_cycle, :r0, 35.0)
    update_param!(m, :co2_cycle, :rC, 0.019)
    update_param!(m, :co2_cycle, :rT, 4.165)
    update_param!(m, :co2_cycle, :iIRF_max, 97.0)
    update_param!(m, :co2_cycle, :a, [0.2173, 0.2240, 0.2824, 0.2763])
    update_param!(m, :co2_cycle, :τ_CO₂, [10.0^6, 394.4, 36.54, 4.304])
    update_param!(m, :co2_cycle, :E_CO₂, rcp_emissions.FossilCO2 .+ rcp_emissions.OtherCO2)

    # ---- Other Well-Mixed Greenhouse Gas Cycles ---- #
    update_param!(m, :other_ghg_cycles, :τ_other_ghg, gas_data[findall((in)(other_ghg_names), gas_data.gas), :lifetimes])    
    update_param!(m, :other_ghg_cycles, :emiss_other_ghg, Matrix(rcp_emissions[!,Symbol.(other_ghg_names)]))
    update_param!(m, :other_ghg_cycles, :emiss2conc_other_ghg, conversions[findall((in)(other_ghg_names), conversions.gases), :emiss2conc])

    # ---- Methane Radiative Forcing ---- #
    update_param!(m, :ch4_rf, :H₂O_share, 0.12)
    update_param!(m, :ch4_rf, :scale_CH₄, 1.0)
    update_param!(m, :ch4_rf, :a₃, -1.3e-6)
    update_param!(m, :ch4_rf, :b₃, -8.2e-6)

    # ---- Nitrous Oxide Radiative Forcing ---- #
    update_param!(m, :n2o_rf, :a₂, -8.0e-6)
    update_param!(m, :n2o_rf, :b₂, 4.2e-6)
    update_param!(m, :n2o_rf, :c₂, -4.9e-6)    

    # ---- Carbon Dioxide Radiative Forcing ---- #
    update_param!(m, :co2_rf, :a₁, -2.4e-7)
    update_param!(m, :co2_rf, :b₁, 7.2e-4)
    update_param!(m, :co2_rf, :c₁, -2.1e-4)
    update_param!(m, :co2_rf, :rf_scale_CO₂, scale_co2_forcing)

    # ---- Other Well-Mixed Greenhouse Gas Radiative Forcings ---- #
    update_param!(m, :other_ghg_rf, :radiative_efficiency, gas_data[findall((in)(other_ghg_names), gas_data.gas), :rad_eff])

    # ---- Tropospheric Ozone Radiative Forcing ---- #
    update_param!(m, :trop_o3_rf, :NOx_emissions, rcp_emissions.NOx)
    update_param!(m, :trop_o3_rf, :CO_emissions, rcp_emissions.CO)
    update_param!(m, :trop_o3_rf, :NMVOC_emissions, rcp_emissions.NMVOC)    
    update_param!(m, :trop_o3_rf, :mol_weight_NO, gas_data[gas_data.gas .== "NO", :mol_weight][1])
    update_param!(m, :trop_o3_rf, :T0, 0.0)    

    # ---- Stratospheric Ozone Radiative Forcing ---- #
    update_param!(m, :strat_o3_rf, :Br, gas_data[findall((in)(ods_names), gas_data.gas), :br_atoms])
    update_param!(m, :strat_o3_rf, :Cl, gas_data[findall((in)(ods_names), gas_data.gas), :cl_atoms])
    update_param!(m, :strat_o3_rf, :FC, gas_data[findall((in)(ods_names), gas_data.gas), :strat_frac])
    update_param!(m, :strat_o3_rf, :δ1, -1.46030698e-5)
    update_param!(m, :strat_o3_rf, :δ2, 2.05401270e-3)
    update_param!(m, :strat_o3_rf, :δ3, 1.03143308)
    update_param!(m, :strat_o3_rf, :ODS₀, gas_data[findall((in)(ods_names), gas_data.gas), :pi_conc])

    # ---- Aerosol Direct Radiative Forcing ---- #
    update_param!(m, :aerosol_direct_rf, :β_SOx, -6.2227e-3)
    update_param!(m, :aerosol_direct_rf, :β_CO, 0.0)
    update_param!(m, :aerosol_direct_rf, :β_NMVOC, -3.8392e-4)
    update_param!(m, :aerosol_direct_rf, :β_NOx, -1.16551e-3)
    update_param!(m, :aerosol_direct_rf, :β_BC, 1.601537e-2)
    update_param!(m, :aerosol_direct_rf, :β_OC, -1.45339e-3)
    update_param!(m, :aerosol_direct_rf, :β_NH3, -1.55605e-3)
    update_param!(m, :aerosol_direct_rf, :rf_scale_aero_direct, 0.0)    
    update_param!(m, :aerosol_direct_rf, :CO_emiss, rcp_emissions.CO)
    update_param!(m, :aerosol_direct_rf, :NMVOC_emiss, rcp_emissions.NMVOC)
    update_param!(m, :aerosol_direct_rf, :NH3_emiss, rcp_emissions.NH3)

    # ---- Aerosol Indirect Radiative Forcing ---- #
    update_param!(m, :aerosol_indirect_rf, :ϕ, -1.95011431)
    update_param!(m, :aerosol_indirect_rf, :b_SOx, 0.01107147)
    update_param!(m, :aerosol_indirect_rf, :b_POM, 0.01387492)
    update_param!(m, :aerosol_indirect_rf, :rf_scale_aero_indirect, 0.0)
    update_param!(m, :aerosol_indirect_rf, :model_years, collect(start_year:end_year))
    update_param!(m, :aerosol_indirect_rf, :SOx_emiss_1765, 1.0)
    update_param!(m, :aerosol_indirect_rf, :BC_OC_emiss_1765, 11.2)
    update_param!(m, :aerosol_indirect_rf, :scale_AR5, true)
    update_param!(m, :aerosol_indirect_rf, :F_1765, -0.3002836449793625)
    update_param!(m, :aerosol_indirect_rf, :F_2011, -1.5236182344467388)

    # ---- Black Carbon on Snow Radiative Forcing ---- #

    # ---- Land Use Change Radiative Forcing ---- #
    update_param!(m, :landuse_rf, :landuse_emiss, rcp_emissions.OtherCO2)

    # ---- Contrails Radiative Forcing ---- #    
    update_param!(m, :contrails_rf, :frac, gas_fractions.nox_aviation)
    update_param!(m, :contrails_rf, :E_ref, 2.946)
    update_param!(m, :contrails_rf, :F_ref, 0.0448)
    update_param!(m, :contrails_rf, :ref_is_NO2, true)
    update_param!(m, :contrails_rf, :mol_weight_NO₂, gas_data[gas_data.gas .== "NO2", :mol_weight][1])    

    # ---- Total Radiative Forcing ---- #
    update_param!(m, :total_rf, :F_volcanic, volcano_forcing)
    update_param!(m, :total_rf, :F_solar, solar_forcing)
    update_param!(m, :total_rf, :F_exogenous, zeros(n_steps))
    update_param!(m, :total_rf, :efficacy_CO₂, 1.0)
    update_param!(m, :total_rf, :efficacy_CH₄, 1.0)
    update_param!(m, :total_rf, :efficacy_CH₄_H₂O, 1.0)
    update_param!(m, :total_rf, :efficacy_N₂O, 1.0)
    update_param!(m, :total_rf, :efficacy_other_ghg, ones(length(other_ghg_names)))
    update_param!(m, :total_rf, :efficacy_trop_O₃, 1.0)
    update_param!(m, :total_rf, :efficacy_strat_O₃, 1.0)
    update_param!(m, :total_rf, :efficacy_aerosol_direct, 1.0)
    update_param!(m, :total_rf, :efficacy_aerosol_indirect, 1.0)
    update_param!(m, :total_rf, :efficacy_bcsnow, 3.0)
    update_param!(m, :total_rf, :efficacy_landuse, 1.0)
    update_param!(m, :total_rf, :efficacy_contrails, 0.0) # Note: Efficacy set to 0.0 to match default settings in Python version of FAIR.

    # ---- Global Temperature Anomaly ---- #
    update_param!(m, :temperature, :d, d)
    update_param!(m, :temperature, :q, q)
    update_param!(m, :temperature, :F2x, F2x)

    # --- Set shared parameters ---- #

    # set_param!(m, :CO₂_0, gas_data[gas_data.gas .== "CO2", :pi_conc][1])
    add_shared_param!(m, :CO₂_0, gas_data[gas_data.gas .== "CO2", :pi_conc][1])
    connect_param!(m, :co2_rf, :CO₂_0, :CO₂_0)
    connect_param!(m, :n2o_rf, :CO₂_0, :CO₂_0)

    # set_param!(m, :CH₄_0, gas_data[gas_data.gas .== "CH4", :pi_conc][1])
    add_shared_param!(m, :CH₄_0, gas_data[gas_data.gas .== "CH4", :pi_conc][1])
    connect_param!(m, :ch4_cycle, :CH₄_0, :CH₄_0)
    connect_param!(m, :ch4_rf, :CH₄_0, :CH₄_0)
    connect_param!(m, :n2o_rf, :CH₄_0, :CH₄_0)
    connect_param!(m, :trop_o3_rf, :CH₄_0, :CH₄_0)

    # set_param!(m, :N₂O_0, gas_data[gas_data.gas .== "N2O", :pi_conc][1])
    add_shared_param!(m, :N₂O_0, gas_data[gas_data.gas .== "N2O", :pi_conc][1])
    connect_param!(m, :ch4_rf, :N₂O_0, :N₂O_0)
    connect_param!(m, :co2_rf, :N₂O_0, :N₂O_0)
    connect_param!(m, :n2o_rf, :N₂O_0, :N₂O_0)
    connect_param!(m, :n2o_cycle, :N₂O_0, :N₂O_0)

    # set_param!(m, :other_ghg_0, gas_data[findall((in)(other_ghg_names), gas_data.gas), :pi_conc])
    add_shared_param!(m, :other_ghg_0, gas_data[findall((in)(other_ghg_names), gas_data.gas), :pi_conc], dims = [:other_ghg])
    connect_param!(m, :other_ghg_cycles, :other_ghg_0, :other_ghg_0)
    connect_param!(m, :other_ghg_rf, :other_ghg_0, :other_ghg_0)

    # set_param!(m, :SOx_emiss, rcp_emissions.SOx)
    add_shared_param!(m, :SOx_emiss, rcp_emissions.SOx, dims = [:time])
    connect_param!(m, :aerosol_direct_rf, :SOx_emiss, :SOx_emiss)
    connect_param!(m, :aerosol_indirect_rf, :SOx_emiss, :SOx_emiss)

    # set_param!(m, :BC_emiss, rcp_emissions.BC)
    add_shared_param!(m, :BC_emiss, rcp_emissions.BC, dims = [:time])
    connect_param!(m, :aerosol_direct_rf, :BC_emiss, :BC_emiss)
    connect_param!(m, :aerosol_indirect_rf, :BC_emiss, :BC_emiss)
    connect_param!(m, :bc_snow_rf, :BC_emiss, :BC_emiss)

    # set_param!(m, :OC_emiss, rcp_emissions.OC)
    add_shared_param!(m, :OC_emiss, rcp_emissions.OC, dims = [:time])
    connect_param!(m, :aerosol_direct_rf, :OC_emiss, :OC_emiss)
    connect_param!(m, :aerosol_indirect_rf, :OC_emiss, :OC_emiss)

    # set_param!(m, :NOx_emiss, rcp_emissions.NOx)
    add_shared_param!(m, :NOx_emiss, rcp_emissions.NOx, dims = [:time])
    connect_param!(m, :aerosol_direct_rf, :NOx_emiss, :NOx_emiss)
    connect_param!(m, :contrails_rf, :NOx_emiss, :NOx_emiss)

    # set_param!(m, :fix_pre1850_RCP, true)
    add_shared_param!(m, :fix_pre1850_RCP, true; data_type = Bool)
    connect_param!(m, :aerosol_indirect_rf, :fix_pre1850_RCP, :fix_pre1850_RCP)
    connect_param!(m, :trop_o3_rf, :fix_pre1850_RCP, :fix_pre1850_RCP)

    # set_param!(m, :mol_weight_N, gas_data[gas_data.gas .== "N", :mol_weight][1])
    add_shared_param!(m, :mol_weight_N, gas_data[gas_data.gas .== "N", :mol_weight][1])
    connect_param!(m, :contrails_rf, :mol_weight_N, :mol_weight_N)
    connect_param!(m, :trop_o3_rf, :mol_weight_N, :mol_weight_N)

    # set_param!(m, :gtc2ppm, conversions[conversions.gases .== "CO2", :emiss2conc][1])
    add_shared_param!(m, :gtc2ppm, conversions[conversions.gases .== "CO2", :emiss2conc][1])
    connect_param!(m, :ch4_cycle, :gtc2ppm, :gtc2ppm)
    connect_param!(m, :co2_cycle, :gtc2ppm, :gtc2ppm)

    # ---------------------------------------------
    # Create connections between Mimi components.
    # ---------------------------------------------
    connect_param!(m, :co2_cycle, :T, :temperature, :T)
    connect_param!(m, :co2_cycle, :E_ox_CH₄, :ch4_cycle, :oxidised_CH₄_GtC)
    connect_param!(m, :trop_o3_rf, :CH₄, :ch4_cycle, :CH₄)
    connect_param!(m, :strat_o3_rf, :conc_ODS, :other_ghg_cycles, :conc_ods)
    connect_param!(m, :ch4_rf, :CH₄, :ch4_cycle, :CH₄)
    connect_param!(m, :ch4_rf, :N₂O, :n2o_cycle, :N₂O)
    connect_param!(m, :n2o_rf, :N₂O, :n2o_cycle, :N₂O)
    connect_param!(m, :n2o_rf, :CH₄, :ch4_cycle, :CH₄)
    connect_param!(m, :n2o_rf, :CO₂, :co2_cycle, :C)
    connect_param!(m, :co2_rf, :N₂O, :n2o_cycle, :N₂O)
    connect_param!(m, :co2_rf, :CO₂, :co2_cycle, :C)
    connect_param!(m, :other_ghg_rf, :conc_other_ghg, :other_ghg_cycles, :conc_other_ghg)
    connect_param!(m, :trop_o3_rf, :temperature, :temperature, :T)
    connect_param!(m, :total_rf, :F_CO₂, :co2_rf, :rf_co2)
    connect_param!(m, :total_rf, :F_CH₄, :ch4_rf, :forcing_CH₄)
    connect_param!(m, :total_rf, :F_CH₄_H₂O, :ch4_rf, :forcing_CH₄_H₂O)
    connect_param!(m, :total_rf, :F_N₂O, :n2o_rf, :forcing_N₂O)
    connect_param!(m, :total_rf, :F_other_ghg, :other_ghg_rf, :other_ghg_rf)
    connect_param!(m, :total_rf, :F_trop_O₃, :trop_o3_rf, :forcing_trop_O₃)
    connect_param!(m, :total_rf, :F_strat_O₃, :strat_o3_rf, :forcing_strat_O₃)
    connect_param!(m, :total_rf, :F_aerosol_direct, :aerosol_direct_rf, :F_aerosol_direct)
    connect_param!(m, :total_rf, :F_aerosol_indirect, :aerosol_indirect_rf, :ERF_aero_cloud)
    connect_param!(m, :total_rf, :F_bcsnow, :bc_snow_rf, :forcing_BC_snow)
    connect_param!(m, :total_rf, :F_landuse, :landuse_rf, :forcing_landuse)
    connect_param!(m, :total_rf, :F_contrails, :contrails_rf, :forcing_contrails)
    connect_param!(m, :temperature, :F, :total_rf, :total_forcing)
end

"""
    update_MimiFAIR162_params!(m; ar6_scenario::String="ssp245", start_year::Int=1750, end_year::Int=2300)

Update all parameter settings in the model `m` with default FAIR  v1.2.6 settings 
as pulled from the MimiFAIRv1_6_2 repository.  Optional arguments include `ar6_scenario` to indicate
what scenario to use, `start_year` and `end_year` to indicate the model time dimension.
"""
function update_MimiFAIR162_params!(m; ar6_scenario::String="ssp245", start_year::Int=1750, end_year::Int=2300)
    
    # ---------------------------------------------
    # Set Up Data and Parameter Values
    # ---------------------------------------------

    # Load RCP and other data needed to construct FAIR (just hard-coded 1765 as start year because not using the RCP emissions from this function).
    rcp_emissions, volcano_forcing, solar_forcing, gas_data, gas_fractions, conversions = MimiFAIRv1_6_2.load_fair_data(1765, end_year, "RCP85")

    # Load IPCC AR6 emissions scenario used for FAIRv1.6.2 ensemble runs (options = "ssp119", "ssp126", "ssp245", "ssp370", "ssp460", "ssp585").
    ar6_emissions_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "AR6_emissions_"*ar6_scenario*"_1750_2300.csv")))

    # Subset AR6 emissions to proper years.
    emission_indices = indexin(collect(start_year:end_year), ar6_emissions_raw.Year)
    ar6_emissions = ar6_emissions_raw[emission_indices, :]

    # Load IPCC AR6 natural CH₄ and N₂O emissions for FAIR (spans 1750-2500).
    ar6_natural_emissions_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "fair_wg3_natural_ch4_n2o.csv")))

    # Subset natural CH₄ and N₂O emissions to proper years.
    natural_indices = indexin(collect(start_year:end_year), ar6_natural_emissions_raw.year)
    ar6_natural_emissions = ar6_natural_emissions_raw[natural_indices, [:CH4, :N2O]]

    # Load IPCC solar forcing scenario (note, dataset runs from -6755 to 2299).
    ar6_solar_forcing_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "ar6_solar_erf.csv")))

    # For cases when you want to run the model out to 2300, assume solar forcing value in 2300 follows the trend from 2298-2299.
    push!(ar6_solar_forcing_raw, [2300 ar6_solar_forcing_raw[end,"solar_erf"] + (ar6_solar_forcing_raw[end,"solar_erf"] - ar6_solar_forcing_raw[(end-1),"solar_erf"])])

    # Subset solar forcing data to proper years.
    solar_indices = indexin(collect(start_year:end_year), ar6_solar_forcing_raw.year)
    ar6_solar_forcing = ar6_solar_forcing_raw[solar_indices, :solar_erf]

    # Load IPCC AR6 volcanic forcing scenario (note, dataset runs from -500 to 2019).
    ar6_volcanic_forcing_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "ar6_volcanic_erf.csv")))

    # Extract indices for relevant years.
    if end_year > 2019
        volcanic_indices = indexin(collect(start_year:2019), ar6_volcanic_forcing_raw.year)
    else
        volcanic_indices = indexin(collect(start_year:end_year), ar6_volcanic_forcing_raw.year)
    end

    # Create an empty array and store subset of volcanic forcings.
    ar6_volcanic_forcing = zeros(length(start_year:end_year))
    ar6_volcanic_forcing[1:length(volcanic_indices)] = ar6_volcanic_forcing_raw[volcanic_indices, :volcanic_erf]

    # From FAIR AR6 code on volcanic forcing: "ramp down last 10 years to zero according to https://www.geosci-model-dev.net/9/3461/2016/gmd-9-3461-2016.html"
    # This copies that code exactly.
    index_2019 = findfirst(x->x==2019, start_year:end_year)
    ar6_volcanic_forcing[index_2019:(index_2019+10)] = ar6_volcanic_forcing[index_2019] .* collect(range(1,0,length=11))

    # Names of minor greenhouse gases and ozone-depleting substances (used or indexing).
    other_ghg_names = ["CF4", "C2F6", "C6F14", "HFC23", "HFC32", "HFC43_10", "HFC125", "HFC134a", "HFC143a", "HFC227ea", "HFC245fa", "SF6"]
    ods_names       = ["CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]

    # ---------------------------------------------
    # Set component-specific parameters
    # ---------------------------------------------

    # ---- Carbon Cycle ---- #
    update_param!(m, :co2_cycle, :CO₂_0,  278.052989189439) # From FAIR model run.
    update_param!(m, :co2_cycle, :iirf_h,  100.0)
    update_param!(m, :co2_cycle, :r0_co2, 35.0)
    update_param!(m, :co2_cycle, :rT_co2, 4.165)
    update_param!(m, :co2_cycle, :rC_co2,  0.019)
    update_param!(m, :co2_cycle, :τ_co2, [1000000, 394.4, 36.54, 4.304])
    update_param!(m, :co2_cycle, :a_co2, [0.2173,0.2240,0.2824,0.2763])
    update_param!(m, :co2_cycle, :R0_co2, [0.0003062168651584551, 0.0003156584344017209, 0.0003979550976564552, 0.0003893590420767655]) # From FAIR model run.
    update_param!(m, :co2_cycle, :E_co2, ar6_emissions.FossilCO2 .+ ar6_emissions.OtherCO2)
    update_param!(m, :co2_cycle, :cumulative_emissions_CO2₀, 0.003)
    update_param!(m, :co2_cycle, :airborne_emissions_CO2₀, 0.0)
    update_param!(m, :co2_cycle, :iIRF_max, 97.0)
    connect_param!(m, :co2_cycle => :temperature, :temperature => :T)

    # ---- Methane Cycle ---- #
    update_param!(m, :ch4_cycle, :fossil_emiss_CH₄, ar6_emissions.CH4)
    update_param!(m, :ch4_cycle, :natural_emiss_CH₄, ar6_natural_emissions.CH4)
    update_param!(m, :ch4_cycle, :τ_CH₄, 9.3)
    update_param!(m, :ch4_cycle, :fossil_frac, ones(length(start_year:end_year)))
    update_param!(m, :ch4_cycle, :oxidation_frac, 0.61)
    update_param!(m, :ch4_cycle, :mol_weight_CH₄, gas_data[gas_data.gas .== "CH4", :mol_weight][1])
    update_param!(m, :ch4_cycle, :mol_weight_C, gas_data[gas_data.gas .== "C", :mol_weight][1])
    update_param!(m, :ch4_cycle, :emiss2conc_ch4, conversions[conversions.gases .== "CH4", :emiss2conc][1])
    update_param!(m, :ch4_cycle, :CH₄_0, gas_data[gas_data.gas .== "CH4", :pi_conc_ar6][1])

    # ---- Nitrous Oxide Cycle ---- #
    update_param!(m, :n2o_cycle, :fossil_emiss_N₂O, ar6_emissions.N2O)
    update_param!(m, :n2o_cycle, :natural_emiss_N₂O, ar6_natural_emissions.N2O)
    update_param!(m, :n2o_cycle, :τ_N₂O, 121.0)
    update_param!(m, :n2o_cycle, :emiss2conc_n2o, conversions[conversions.gases .== "N2O", :emiss2conc][1])
    update_param!(m, :n2o_cycle, :N₂O_0, gas_data[gas_data.gas .== "N2O", :pi_conc_ar6][1])

    # ---- Other Well-Mixed Greenhouse Gas Cycles ---- #
    update_param!(m, :other_ghg_cycles, :τ_other_ghg, gas_data[findall((in)(other_ghg_names), gas_data.gas), :lifetimes])
    update_param!(m, :other_ghg_cycles, :emiss_other_ghg, Matrix(ar6_emissions[!,Symbol.(other_ghg_names)]))
    update_param!(m, :other_ghg_cycles, :emiss2conc_other_ghg, conversions[findall((in)(other_ghg_names), conversions.gases), :emiss2conc])
    update_param!(m, :other_ghg_cycles, :other_ghg_0, gas_data[findall((in)(other_ghg_names), gas_data.gas), :pi_conc_ar6])

    # ---- Ozone-Depleting Substance Gas Cycles ---- #
    update_param!(m, :o3_depleting_substance_cycles, :τ_ods, gas_data[findall((in)(ods_names), gas_data.gas), :lifetimes])
    update_param!(m, :o3_depleting_substance_cycles, :emiss_ods, Matrix(ar6_emissions[!,Symbol.(ods_names)]))
    update_param!(m, :o3_depleting_substance_cycles, :emiss2conc_ods, conversions[findall((in)(ods_names), conversions.gases), :emiss2conc])
    update_param!(m, :o3_depleting_substance_cycles, :ods_0, gas_data[findall((in)(ods_names), gas_data.gas), :pi_conc_ar6])

    # ---- Carbon Dioxide Radiative Forcing ---- #
    update_param!(m, :co2_forcing, :F2x, 3.71)
    update_param!(m, :co2_forcing, :a₁, -2.4785e-07)
    update_param!(m, :co2_forcing, :b₁, 0.00075906)
    update_param!(m, :co2_forcing, :c₁, -0.0021492)
    update_param!(m, :co2_forcing, :d₁, 5.2488)
    update_param!(m, :co2_forcing, :adjust_F2x, true)

    connect_param!(m, :co2_forcing => :CO₂, :co2_cycle => :co2)
    connect_param!(m, :co2_forcing => :N₂O, :n2o_cycle => :N₂O)

    # ---- Methane Radiative Forcing ---- #
    update_param!(m, :ch4_forcing, :a₃, -8.9603e-05)
    update_param!(m, :ch4_forcing, :b₃, -0.00012462)
    update_param!(m, :ch4_forcing, :d₃, 0.045194)
    update_param!(m, :ch4_forcing, :h2o_from_ch4, 0.079047)

    connect_param!(m, :ch4_forcing => :N₂O, :n2o_cycle => :N₂O)
    connect_param!(m, :ch4_forcing => :CH₄, :ch4_cycle => :CH₄)

    # ---- Nitrous Oxide Radiative Forcing ---- #
    update_param!(m, :n2o_forcing, :a₂, -0.00034197)
    update_param!(m, :n2o_forcing, :b₂,  0.00025455)
    update_param!(m, :n2o_forcing, :c₂, -0.00024357)
    update_param!(m, :n2o_forcing, :d₂, 0.12173)

    connect_param!(m, :n2o_forcing => :CO₂, :co2_cycle => :co2)
    connect_param!(m, :n2o_forcing => :N₂O, :n2o_cycle => :N₂O)
    connect_param!(m, :n2o_forcing => :CH₄, :ch4_cycle => :CH₄)

    # ---- Ozone Radiative Forcing ---- #
    update_param!(m, :o3_forcing, :total_forcing_O₃_0, 0.0)
    update_param!(m, :o3_forcing, :Br, gas_data[findall((in)(ods_names), gas_data.gas), :br_atoms])
    update_param!(m, :o3_forcing, :Cl, gas_data[findall((in)(ods_names), gas_data.gas), :cl_atoms])
    update_param!(m, :o3_forcing, :FC, gas_data[findall((in)(ods_names), gas_data.gas), :strat_frac])
    update_param!(m, :o3_forcing, :feedback, -0.037)
    update_param!(m, :o3_forcing, :Ψ_CH₄, 2.33379720e-04)
    update_param!(m, :o3_forcing, :Ψ_N₂O, 1.27179106e-03)
    update_param!(m, :o3_forcing, :Ψ_ODS, -6.69347820e-05)
    update_param!(m, :o3_forcing, :Ψ_CO, 1.14647701e-04)
    update_param!(m, :o3_forcing, :Ψ_NMVOC, 5.14366051e-12)
    update_param!(m, :o3_forcing, :Ψ_NOx, 3.78354423e-03)
    update_param!(m, :o3_forcing, :CO_emiss_pi,    348.527359)
    update_param!(m, :o3_forcing, :NMVOC_emiss_pi, 60.0218262)
    update_param!(m, :o3_forcing, :NOx_emiss_pi,   3.87593407)
    
    connect_param!(m, :o3_forcing => :N₂O, :n2o_cycle => :N₂O)
    connect_param!(m, :o3_forcing => :CH₄, :ch4_cycle => :CH₄)
    connect_param!(m, :o3_forcing => :temperature, :temperature => :T)
    connect_param!(m, :o3_forcing => :conc_ODS, :o3_depleting_substance_cycles => :conc_ods)

    # ---- Aerosol Direct Radiative Forcing ---- #
    update_param!(m, :aerosol_direct_forcing, :β_SOx, -6.2227e-3)
    update_param!(m, :aerosol_direct_forcing, :β_CO, 0.0)
    update_param!(m, :aerosol_direct_forcing, :β_NMVOC, -3.8392e-4)
    update_param!(m, :aerosol_direct_forcing, :β_NOx, -1.16551e-3)
    update_param!(m, :aerosol_direct_forcing, :β_BC, 1.601537e-2)
    update_param!(m, :aerosol_direct_forcing, :β_OC, -1.45339e-3)
    update_param!(m, :aerosol_direct_forcing, :β_NH3, -1.55605e-3)
    update_param!(m, :aerosol_direct_forcing, :NH3_emiss, ar6_emissions.NH3)

    # ---- Aerosol Indirect Radiative Forcing ---- #
    update_param!(m, :aerosol_indirect_forcing, :ϕ, 0.07334277994353743)#-1.95011431)
    update_param!(m, :aerosol_indirect_forcing, :b_SOx, 3.452849302362568)#0.01107147)
    update_param!(m, :aerosol_indirect_forcing, :b_POM, 33.126485122209154)#0.01387492)
    update_param!(m, :aerosol_indirect_forcing, :rf_scale_aero_indirect, 1.0)
    update_param!(m, :aerosol_indirect_forcing, :SOx_emiss_pi,   1.22002422)
    update_param!(m, :aerosol_indirect_forcing, :BC_emiss_pi,    2.09777075)
    update_param!(m, :aerosol_indirect_forcing, :OC_emiss_pi,    15.4476682)

    # ---- Other Well-Mixed Greenhouse Gas Radiative Forcings ---- #
    update_param!(m, :other_ghg_forcing, :other_ghg_radiative_efficiency, gas_data[findall((in)(other_ghg_names), gas_data.gas), :rad_eff])
    update_param!(m, :other_ghg_forcing, :other_ghg_pi, gas_data[findall((in)(other_ghg_names), gas_data.gas), :pi_conc_ar6])

    connect_param!(m, :other_ghg_forcing => :conc_other_ghg, :other_ghg_cycles => :conc_other_ghg)

    # ---- Ozone-Depleting Substance Radiative Forcings ---- #
    update_param!(m, :o3_depleting_substance_forcing, :ods_radiative_efficiency, gas_data[findall((in)(ods_names), gas_data.gas), :rad_eff])

    connect_param!(m, :o3_depleting_substance_forcing => :conc_ods, :o3_depleting_substance_cycles => :conc_ods)

    # ---- Contrails Radiative Forcing ---- #
    update_param!(m, :contrails_forcing, :frac, zeros(length(start_year:end_year)))
    update_param!(m, :contrails_forcing, :E_ref_contrails, 2.946)
    update_param!(m, :contrails_forcing, :F_ref_contrails, 0.0448)
    update_param!(m, :contrails_forcing, :ref_is_NO2, true)
    update_param!(m, :contrails_forcing, :mol_weight_NO₂, gas_data[gas_data.gas .== "NO2", :mol_weight][1])
    update_param!(m, :contrails_forcing, :mol_weight_N, gas_data[gas_data.gas .== "N", :mol_weight][1])

    # ---- Black Carbon on Snow Radiative Forcing ---- #
    update_param!(m, :bc_snow_forcing, :E_ref_bc, 6.095)
    update_param!(m, :bc_snow_forcing, :F_ref_bc, 0.08)

    # ---- Land Use Change Radiative Forcing ---- #
    update_param!(m, :landuse_forcing, :α_CO₂_land, (-0.2/190))
    update_param!(m, :landuse_forcing, :landuse_emiss, ar6_emissions.OtherCO2)

    # ---- Total Radiative Forcing ---- #
    update_param!(m, :total_forcing, :scale_CO₂, 1.0)
    update_param!(m, :total_forcing, :scale_CH₄, 1.0)
    update_param!(m, :total_forcing, :scale_CH₄_H₂O, 1.0)
    update_param!(m, :total_forcing, :scale_N₂O, 1.0)
    update_param!(m, :total_forcing, :scale_O₃, 1.0)
    update_param!(m, :total_forcing, :scale_aerosol_indirect, 1.0)
    update_param!(m, :total_forcing, :scale_bcsnow, 1.0)
    update_param!(m, :total_forcing, :scale_landuse, 1.0)
    update_param!(m, :total_forcing, :scale_contrails, 0.0) # Default FAIR has contrail forcing switched off. Set scaling term to 0
    update_param!(m, :total_forcing, :scale_volcanic, 1.0)
    update_param!(m, :total_forcing, :scale_solar, 1.0)
    update_param!(m, :total_forcing, :scale_aerosol_direct_SOx, 1.0)
    update_param!(m, :total_forcing, :scale_aerosol_direct_CO_NMVOC, 1.0)
    update_param!(m, :total_forcing, :scale_aerosol_direct_NOx_NH3, 1.0)
    update_param!(m, :total_forcing, :scale_aerosol_direct_BC, 1.0)
    update_param!(m, :total_forcing, :scale_aerosol_direct_OC, 1.0)
    update_param!(m, :total_forcing, :scale_other_ghg, ones(length(other_ghg_names)))
    update_param!(m, :total_forcing, :scale_ods, ones(length(ods_names)))
    update_param!(m, :total_forcing, :F_volcanic, ar6_volcanic_forcing)
    update_param!(m, :total_forcing, :F_solar, ar6_solar_forcing)
    update_param!(m, :total_forcing, :F_exogenous, zeros(length(start_year:end_year)))

    connect_param!(m, :total_forcing => :F_CO₂, :co2_forcing => :rf_co2)
    connect_param!(m, :total_forcing => :F_CH₄, :ch4_forcing => :rf_ch4)
    connect_param!(m, :total_forcing => :F_CH₄_H₂O, :ch4_forcing => :rf_ch4_h2o)
    connect_param!(m, :total_forcing => :F_N₂O, :n2o_forcing => :rf_n2o)
    connect_param!(m, :total_forcing => :F_O₃, :o3_forcing => :total_forcing_O₃)
    connect_param!(m, :total_forcing => :F_aerosol_direct_SOx, :aerosol_direct_forcing => :F_SOx_aero)
    connect_param!(m, :total_forcing => :F_aerosol_direct_CO, :aerosol_direct_forcing => :F_CO_aero)
    connect_param!(m, :total_forcing => :F_aerosol_direct_NMVOC, :aerosol_direct_forcing => :F_NMVOC_aero)
    connect_param!(m, :total_forcing => :F_aerosol_direct_NOx, :aerosol_direct_forcing => :F_NOx_aero)
    connect_param!(m, :total_forcing => :F_aerosol_direct_BC, :aerosol_direct_forcing => :F_BC_aero)
    connect_param!(m, :total_forcing => :F_aerosol_direct_OC, :aerosol_direct_forcing => :F_OC_aero)
    connect_param!(m, :total_forcing => :F_aerosol_direct_NH3, :aerosol_direct_forcing => :F_NH3_aero)
    connect_param!(m, :total_forcing => :F_aerosol_indirect, :aerosol_indirect_forcing => :rf_aero_indirect)
    connect_param!(m, :total_forcing => :F_bcsnow, :bc_snow_forcing => :forcing_BC_snow)
    connect_param!(m, :total_forcing => :F_landuse, :landuse_forcing => :forcing_landuse)
    connect_param!(m, :total_forcing => :F_contrails, :contrails_forcing => :forcing_contrails)
    connect_param!(m, :total_forcing => :F_other_ghg, :other_ghg_forcing => :other_ghg_rf)
    connect_param!(m, :total_forcing => :F_ods, :o3_depleting_substance_forcing => :ods_rf)

    # ---- Temperature ---- #
    update_param!(m, :temperature, :earth_radius, 6371000)
    update_param!(m, :temperature, :seconds_per_year, (60*60*24*365.24219))
    update_param!(m, :temperature, :ocean_heat_exchange, 0.67)
    update_param!(m, :temperature, :deep_ocean_efficacy, 1.28)
    update_param!(m, :temperature, :lambda_global, 1.18)
    update_param!(m, :temperature, :T_mix₀, [5.30299681e-03, 6.41290105e-05]) # From Python-FAIR model run.
    update_param!(m, :temperature, :T_deep₀, [-1.33065374e-04,  1.50206328e-04]) # From Python-FAIR model run.
    update_param!(m, :temperature, :ocean_heat_capacity, [8.2, 109.0])
    connect_param!(m, :temperature => :forcing, :total_forcing => :total_forcing)

    # ---- Parameters Shared Across Multiple Components - Set shared parameters ---- #
    
    # set_param!(m, :dt, 1.0)
    add_shared_param!(m, :model_dt, 1.0)
    connect_param!(m, :co2_cycle,   :dt, :model_dt)
    connect_param!(m, :temperature, :dt, :model_dt)

    # set_param!(m, :emiss2conc_co2, conversions[conversions.gases .== "CO2", :emiss2conc][1])
    add_shared_param!(m, :model_emiss2conc_co2, conversions[conversions.gases .== "CO2", :emiss2conc][1])
    connect_param!(m, :co2_cycle, :emiss2conc_co2, :model_emiss2conc_co2)
    connect_param!(m, :ch4_cycle, :emiss2conc_co2, :model_emiss2conc_co2)

    # set_param!(m, :CH₄_pi, gas_data[gas_data.gas .== "CH4", :pi_conc_ar6][1])
    add_shared_param!(m, :model_CH₄_pi, gas_data[gas_data.gas .== "CH4", :pi_conc_ar6][1])
    connect_param!(m, :ch4_cycle,   :CH₄_pi, :model_CH₄_pi)
    connect_param!(m, :ch4_forcing, :CH₄_pi, :model_CH₄_pi)
    connect_param!(m, :o3_forcing,  :CH₄_pi, :model_CH₄_pi)

    # set_param!(m, :CO₂_pi, gas_data[gas_data.gas .== "CO2", :pi_conc_ar6][1])
    add_shared_param!(m, :model_CO₂_pi, gas_data[gas_data.gas .== "CO2", :pi_conc_ar6][1])
    connect_param!(m, :co2_cycle,   :CO₂_pi, :model_CO₂_pi)
    connect_param!(m, :co2_forcing, :CO₂_pi, :model_CO₂_pi)

    # set_param!(m, :N₂O_pi, gas_data[gas_data.gas .== "N2O", :pi_conc_ar6][1])
    add_shared_param!(m, :model_N₂O_pi, gas_data[gas_data.gas .== "N2O", :pi_conc_ar6][1])
    connect_param!(m, :co2_forcing, :N₂O_pi, :model_N₂O_pi)
    connect_param!(m, :n2o_forcing, :N₂O_pi, :model_N₂O_pi)
    connect_param!(m, :o3_forcing,  :N₂O_pi, :model_N₂O_pi)

    # set_param!(m, :ods_pi, gas_data[findall((in)(ods_names), gas_data.gas), :pi_conc_ar6])
    add_shared_param!(m, :model_ods_pi, gas_data[findall((in)(ods_names), gas_data.gas), :pi_conc_ar6], dims = [:ozone_depleting_substances])
    connect_param!(m, :o3_depleting_substance_forcing,  :ods_pi, :model_ods_pi)
    connect_param!(m, :o3_forcing,                      :ods_pi, :model_ods_pi)

    # set_param!(m, :SOx_emiss, ar6_emissions.SOx)
    add_shared_param!(m, :model_SOx_emiss, ar6_emissions.SOx, dims = [:time])
    connect_param!(m, :aerosol_direct_forcing,      :SOx_emiss, :model_SOx_emiss)
    connect_param!(m, :aerosol_indirect_forcing,    :SOx_emiss, :model_SOx_emiss)

    # set_param!(m, :BC_emiss, ar6_emissions.BC)
    add_shared_param!(m, :model_BC_emiss, ar6_emissions.BC, dims = [:time])
    connect_param!(m, :aerosol_direct_forcing,      :BC_emiss, :model_BC_emiss)
    connect_param!(m, :aerosol_indirect_forcing,    :BC_emiss, :model_BC_emiss)
    connect_param!(m, :bc_snow_forcing,             :BC_emiss, :model_BC_emiss)

    # set_param!(m, :OC_emiss, ar6_emissions.OC)
    add_shared_param!(m, :model_OC_emiss, ar6_emissions.OC, dims = [:time])
    connect_param!(m, :aerosol_direct_forcing,      :OC_emiss, :model_OC_emiss)
    connect_param!(m, :aerosol_indirect_forcing,    :OC_emiss, :model_OC_emiss)

    # set_param!(m, :CO_emiss, ar6_emissions.CO)
    add_shared_param!(m, :model_CO_emiss, ar6_emissions.CO, dims = [:time])
    connect_param!(m, :aerosol_direct_forcing,  :CO_emiss, :model_CO_emiss)
    connect_param!(m, :o3_forcing,              :CO_emiss, :model_CO_emiss)

    # set_param!(m, :NMVOC_emiss, ar6_emissions.NMVOC)
    add_shared_param!(m, :model_NMVOC_emiss, ar6_emissions.NMVOC, dims = [:time])
    connect_param!(m, :aerosol_direct_forcing,  :NMVOC_emiss, :model_NMVOC_emiss)
    connect_param!(m, :o3_forcing,              :NMVOC_emiss, :model_NMVOC_emiss)

    # set_param!(m, :NOx_emiss, ar6_emissions.NOx)
    add_shared_param!(m, :model_NOx_emiss, ar6_emissions.NOx, dims = [:time])
    connect_param!(m, :aerosol_direct_forcing,  :NOx_emiss, :model_NOx_emiss)
    connect_param!(m, :o3_forcing,              :NOx_emiss, :model_NOx_emiss)
    connect_param!(m, :contrails_forcing,       :NOx_emiss, :model_NOx_emiss)

end
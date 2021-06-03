using Mimi
using DataFrames

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
    update_param!(m, :other_ghg_cycles, :emiss_other_ghg, convert(Matrix, rcp_emissions[!,Symbol.(other_ghg_names)]))
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

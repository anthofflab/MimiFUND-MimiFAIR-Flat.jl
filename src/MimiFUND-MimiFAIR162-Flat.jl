using MimiFAIRv1_6_2 # imports the model and its components, as referenced by MimiFAIRv1_6_2.component_name
using MimiFUND
using Mimi
using DataFrames
using VegaLite

include(joinpath(@__DIR__, "helper.jl")) # path relative to this file

function get_fundfair162(;ar6_scenario::String = "ssp245",
                            FAIR_first = 1750,
                            FAIR_last = 2300,
                            FUND_first = 1950,
                            FUND_last = 2300)

    FAIR_len = length(FAIR_first:FAIR_last)
    FUND_len = length(FUND_first:FUND_last)

    # start with default MimiFUND model
    m = MimiFUND.get_model()

    # add new dimensions relevant to FAIR (1) minor greenhouse gases and 
    # (2) ozone-depleting substances
    other_ghg_names = ["CF4", "C2F6", "C6F14", "HFC23", "HFC32", "HFC43_10", "HFC125", "HFC134a", "HFC143a", "HFC227ea", "HFC245fa", "SF6"]
    ods_names       = ["CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]

    set_dimension!(m, :other_ghg, other_ghg_names) # FAIR dims - Set index for Kyoto and ozone ozone-depleting gases.
    set_dimension!(m, :ozone_depleting_substances, ods_names) # FAIR dims - Set index for Kyoto and ozone ozone-depleting gases.

    # set the time dimension of the model to match the full extent of FAIR 
    set_dimension!(m, :time, collect(FAIR_first:FAIR_last))

    # add the MimiFAIR components to the MimiFUND model after the emissions 
    # component, noting they will take the first and last of the model which is now
    # FAIR_first and FAIR_last

    add_comp!(m, MimiFAIRv1_6_2.ch4_cycle; after = :emissions);
    add_comp!(m, MimiFAIRv1_6_2.n2o_cycle; after = :ch4_cycle);
    add_comp!(m, MimiFAIRv1_6_2.co2_cycle; after = :n2o_cycle);
    add_comp!(m, MimiFAIRv1_6_2.other_ghg_cycles; after = :co2_cycle);
    add_comp!(m, MimiFAIRv1_6_2.o3_depleting_substance_cycles; after = :other_ghg_cycles);
    add_comp!(m, MimiFAIRv1_6_2.co2_forcing; after = :o3_depleting_substance_cycles);
    add_comp!(m, MimiFAIRv1_6_2.ch4_forcing; after = :co2_forcing);
    add_comp!(m, MimiFAIRv1_6_2.n2o_forcing; after = :ch4_forcing);
    add_comp!(m, MimiFAIRv1_6_2.o3_forcing; after = :n2o_forcing);
    add_comp!(m, MimiFAIRv1_6_2.aerosol_direct_forcing; after = :o3_forcing);
    add_comp!(m, MimiFAIRv1_6_2.aerosol_indirect_forcing; after = :aerosol_direct_forcing);
    add_comp!(m, MimiFAIRv1_6_2.other_ghg_forcing; after = :aerosol_indirect_forcing);
    add_comp!(m, MimiFAIRv1_6_2.o3_depleting_substance_forcing; after = :other_ghg_forcing);
    add_comp!(m, MimiFAIRv1_6_2.contrails_forcing; after = :o3_depleting_substance_forcing);
    add_comp!(m, MimiFAIRv1_6_2.bc_snow_forcing; after = :contrails_forcing);
    add_comp!(m, MimiFAIRv1_6_2.landuse_forcing; after = :bc_snow_forcing);
    add_comp!(m, MimiFAIRv1_6_2.total_forcing; after = :landuse_forcing);
    add_comp!(m, MimiFAIRv1_6_2.temperature; after = :total_forcing);

    # set all the FAIR component parameters and make their internal connections
    update_MimiFAIR162_params!(m; start_year = FAIR_first, end_year = FAIR_last, ar6_scenario = ar6_scenario)

    # connect FUND and FAIR

    #
    # FUND Emissions --> Convert Mtons to Gtons --> FAIR CO₂ Cycle 
    #

    # destination: FAIR component :co2_cycle parameter :E_CO₂
    # previous source: exogenous parameter
    # new source: FUND :emissions variable :mco2 (which first runs through a multiplier component)

    add_comp!(m, Mimi.multiplier; after = :emissions, first = FUND_first, last = FUND_last);
    set_param!(m, :multiplier, :multiply, fill(1/1000, FAIR_len)) # convert Mtons C coming out of FUND to Gtons C going into FAIR
    connect_param!(m, :multiplier, :input, :emissions, :mco2)

    # Subset AR6 emissions to proper years.
    ar6_emissions_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "AR6_emissions_"*ar6_scenario*"_1750_2300.csv")))
    # ar6_emissions_raw = DataFrame(load(joinpath(@__DIR__, "data", "model_data", "AR6_emissions_"*ar6_scenario*"_1750_2300.csv")))
    emission_indices = indexin(collect(FAIR_first:FAIR_last), ar6_emissions_raw.Year)
    ar6_emissions = ar6_emissions_raw[emission_indices, :]

    FAIR_CO₂_backup = (ar6_emissions.FossilCO2 .+ ar6_emissions.OtherCO2)
    connect_param!(m, :co2_cycle, :E_co2, :multiplier, :output, FAIR_CO₂_backup, backup_offset = 1)

    # A note here is that FAIR v1.6.2 splits land use CO2 emissions and fossil CO2 emission, and feeds 
    # their sum into the `co2_cycle`, but just the land use emissions into the `landuse_forcing` 
    # component in MimiFAIRv1_6_2. Given FUND does not explicilty provide the land use emissions 
    # (AFOLU), we do not modify the latter connection, and thus FAIR will run with its default 
    # land use emissions for the `landuse_forcing` component.

    #
    # FAIR Temperature --> FUND CO₂ Cycle, Biodiversity, Ocean and Regional Climate Components
    #

    # destination: FUND components :climateco2cycle, :biodiversity, :ocean --> parameter :temp
    # previous source: FUND component :climatedynamics --> variable :temp
    # new source: FAIR :temperature --> variable :T
    connect_param!(m, :climateco2cycle, :temp, :temperature, :T)
    connect_param!(m, :biodiversity, :temp, :temperature, :T)
    connect_param!(m, :ocean, :temp, :temperature, :T)

    # destination: FUND component :climateregional, parameter :inputtemp
    # previous source: FUND component :climatedynamics --> variable :temp
    # new source: FAIR component :temperature --> variable :T
    connect_param!(m, :climateregional, :inputtemp, :temperature, :T)

    return m

end
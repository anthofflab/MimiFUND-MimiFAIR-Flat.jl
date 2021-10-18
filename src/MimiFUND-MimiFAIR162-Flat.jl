using MimiFAIRv1_6_2
using MimiFUND
using Mimi

include(joinpath(@__DIR__, "helper.jl")) # path relative to this file

# set some constants
ar6_scenario = "ssp245"

FAIR_first = 1750
FAIR_last = 2300
FAIR_len = length(FAIR_first:FAIR_last)

FUND_first = 1950
FUND_last = 2300
FUND_len = length(FUND_first:FUND_last)

# start with default MimiFUND model
m = MimiFUND.get_model()

# import MimiFAIRv1_6_2 components 
import MimiFAIRv1_6_2: ch4_cycle, n2o_cycle, other_ghg_cycles, co2_cycle, ch4_forcing, o3_forcing,
    n2o_forcing, other_ghg_forcing, other_ghg_cycles, co2_forcing, o3_depleting_substance_cycles, o3_depleting_substance_forcing, aerosol_direct_forcing, 
    aerosol_indirect_forcing, bc_snow_forcing, landuse_forcing, contrails_forcing, total_forcing, temperature 

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
add_comp!(m, ch4_cycle; after = :emissions);
add_comp!(m, n2o_cycle; after = :ch4_cycle);
add_comp!(m, co2_cycle; after = :n2o_cycle);
add_comp!(m, other_ghg_cycles; after = :co2_cycle);
add_comp!(m, o3_depleting_substance_cycles; after = :other_ghg_cycles);
add_comp!(m, co2_forcing; after = :o3_depleting_substance_cycles);
add_comp!(m, ch4_forcing; after = :co2_cycle);
add_comp!(m, n2o_forcing; after = :ch4_forcing);
add_comp!(m, o3_forcing; after = :n2o_forcing);
add_comp!(m, aerosol_direct_forcing; after = :o3_forcing);
add_comp!(m, aerosol_indirect_forcing; after = :aerosol_direct_forcing);
add_comp!(m, other_ghg_forcing; after = :aerosol_indirect_forcing);
add_comp!(m, o3_depleting_substance_forcing; after = :other_ghg_forcing);
add_comp!(m, contrails_forcing; after = :o3_depleting_substance_forcing);
add_comp!(m, bc_snow_forcing; after = :contrails_forcing);
add_comp!(m, landuse_forcing; after = :bc_snow_forcing);
add_comp!(m, total_forcing; after = :landuse_forcing);
add_comp!(m, temperature; after = :total_forcing);

# set all the FAIR component parameters and make their internal connections
set_MimiFAIR_params!(m)

# connect FUND and FAIR

#
# FUND Emissions --> Convert Mtons to Gtons --> FAIR CO₂ Cycle 
#

# destination: FAIR component :co2_cycle parameter :E_CO₂
# previous source: exogenous parameter
# new source: FUND :emissions variable :mco2 (which first runs through a multiplier component)

add_comp!(m, Mimi.multiplier; after = :emissions, first = FUND_first, last = FUND_last);
set_param!(m, :multiplier, :multiply, fill(1/1000, FAIR_len)) # convert Mtons CO₂ coming out of FUND to Gtons CO₂ going into FAIR -- SHOULD THIS BE GTC?
connect_param!(m, :multiplier, :input, :emissions, :mco2)

rcp_emissions, volcano_forcing, solar_forcing, gas_data, gas_fractions, conversions = MimiFAIRv1_6_2.load_fair_data(1765, FAIR_last, "RCP85");

ar6_emissions_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "AR6_emissions_"*ar6_scenario*"_1750_2300.csv")))
# ar6_emissions_raw = DataFrame(load(joinpath(@__DIR__, "data", "model_data", "AR6_emissions_"*ar6_scenario*"_1750_2300.csv")))
# Subset AR6 emissions to proper years.
emission_indices = indexin(collect(FAIR_first:FAIR_last), ar6_emissions_raw.Year)
ar6_emissions = ar6_emissions_raw[emission_indices, :]

FAIR_CO₂_backup = (ar6_emissions.FossilCO2 .+ ar6_emissions.OtherCO2)
connect_param!(m, :co2_cycle, :E_co2, :multiplier, :output, FAIR_CO₂_backup, backup_offset = 1)

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

##
## Step 5. Run and Explore
##

run(m)
explore(m)

fairvals = mfair[:temperature, :T]
fundfairvals =  m[:temperature, :T]
fundvals = vcat(
    fill(missing, length(FAIR_first:FUND_first)-1),
    mfund[:climateco2cycle, :temp][1:length(FUND_first:FAIR_last),:]
)
fundvals = fundvals[:,1]; # make a vector

df = DataFrame(
    :Year => Mimi.time_labels(m),
    :FUND => fundvals,
    :FAIR => fairvals,
    :FUNDFAIR => fundfairvals
)

stack(df, [:FUND, :FAIR, :FUNDFAIR]) |> 

@vlplot(
    :line, 
    x = :Year,
    y = :value,
    color = :variable
)
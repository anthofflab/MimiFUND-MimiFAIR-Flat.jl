using MimiFAIR
using MimiFUND
using Mimi

include("helper.jl")
# include("src/helper.jl")

# set some constants
rcp = "RCP85"

FAIR_first = 1765
FAIR_last = 2500
FAIR_len = length(FAIR_first:FAIR_last)

FUND_first = 1950
FUND_last = 2500
FUND_len = length(FUND_first:FUND_last)

# start with default MimiFUND model
m = MimiFUND.get_model()

# import MimiFAIR components 
import MimiFAIR: ch4_cycle, n2o_cycle, other_ghg_cycles, co2_cycle, ch4_rf, 
    n2o_rf, other_ghg_rf, co2_rf, trop_o3_rf, strat_o3_rf, aerosol_direct_rf, 
    aerosol_indirect_rf, bc_snow_rf, landuse_rf, contrails_rf, total_rf, temperature 

# add new dimensions relevant to FAIR (1) minor greenhouse gases and 
# (2) ozone-depleting substances
other_ghg_names = ["CF4", "C2F6", "C6F14", "HFC23", "HFC32", "HFC43_10", "HFC125", "HFC134a", "HFC143a", "HFC227ea", "HFC245fa", "SF6", "CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]
ods_names = ["CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]

set_dimension!(m, :other_ghg, other_ghg_names) # FAIR dims - Set index for Kyoto and ozone ozone-depleting gases.
set_dimension!(m, :ozone_depleting_substances, ods_names) # FAIR dims - Set index for Kyoto and ozone ozone-depleting gases.

# expand the time dimension of the model to match the full extent of FAIR and 
set_dimension!(m, :time, collect(FAIR_first:FAIR_last))

# add the MimiFAIR components to the MimiFUND model after the emissions 
# component, noting they will take the first and last of the model which is now
# FAIR_first and FAIR_last
add_comp!(m, ch4_cycle; after = :emissions);
add_comp!(m, n2o_cycle; after = :ch4_cycle);
add_comp!(m, other_ghg_cycles; after = :n2o_cycle);
add_comp!(m, co2_cycle; after = :other_ghg_cycles);
add_comp!(m, ch4_rf; after = :co2_cycle);
add_comp!(m, n2o_rf; after = :ch4_rf);
add_comp!(m, other_ghg_rf; after = :n2o_rf);
add_comp!(m, co2_rf; after = :other_ghg_rf);
add_comp!(m, trop_o3_rf; after = :co2_rf);
add_comp!(m, strat_o3_rf; after = :trop_o3_rf);
add_comp!(m, aerosol_direct_rf; after = :strat_o3_rf);
add_comp!(m, aerosol_indirect_rf; after = :aerosol_direct_rf);
add_comp!(m, bc_snow_rf; after = :aerosol_indirect_rf);
add_comp!(m, landuse_rf; after = :bc_snow_rf);
add_comp!(m, contrails_rf; after = :landuse_rf);
add_comp!(m, total_rf; after = :contrails_rf);
add_comp!(m, temperature; after = :total_rf);

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
set_param!(m, :multiplier, :multiply, fill(1/1000, FAIR_len)) # convert Mtons CO₂ coming out of FUND to Gtons CO₂ going into FAIR
connect_param!(m, :multiplier, :input, :emissions, :mco2)

rcp_emissions, volcano_forcing, solar_forcing, gas_data, gas_fractions, conversions = MimiFAIR.load_fair_data(FAIR_first, FAIR_last, rcp);
FAIR_CO₂_backup = (rcp_emissions.FossilCO2 .+ rcp_emissions.OtherCO2)
connect_param!(m, :co2_cycle, :E_CO₂, :multiplier, :output, FAIR_CO₂_backup, backup_offset = 1)

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
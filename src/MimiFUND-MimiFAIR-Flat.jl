using MimiFAIR
using MimiFUND
using Mimi

include("helper.jl")

# set some constants
FAIR_first = 1765
FAIR_last = 2500
FUND_first = 1950
FUND_last = 2500

rcp = "RCP85"

# start with default MimiFUND model and grab a list of the component names, which 
# was originally generated with [keys(m.md.namespace)...] on the FUND model 
m = MimiFUND.get_model()
FUND_comp_names = [:scenariouncertainty, :population,:geography,:socioeconomic,
    :emissions,:climateco2cycle, :climatech4cycle,:climaten2ocycle,:climatesf6cycle,
    :climateforcing,:climatedynamics, :biodiversity,:climateregional,:ocean,
    :impactagriculture,:impactbiodiversity,:impactcardiovascularrespiratory,
    :impactcooling,:impactdiarrhoea,:impactextratropicalstorms,:impactforests,
    :impactheating,:impactvectorbornediseases,:impacttropicalstorms,:vslvmorb,
    :impactdeathmorbidity,:impactwaterresources,:impactsealevelrise,:impactaggregation]

# import MimiFAIR components and our own co2_cycle component with the units 
# adjustment so it properly matches FUND
import MimiFAIR: ch4_cycle, n2o_cycle, other_ghg_cycles, ch4_rf, 
    n2o_rf, other_ghg_rf, co2_rf, trop_o3_rf, strat_o3_rf, aerosol_direct_rf, 
    aerosol_indirect_rf, bc_snow_rf, landuse_rf, contrails_rf, total_rf, temperature 

include("MimiFAIR_co2_cycle_units.jl") # TODO unit conversion component

# add new dimensions relevant to FAIR (1) regions (2) minor greenhouse gases and 
# (3) ozone-depleting substances
other_ghg_names = ["CF4", "C2F6", "C6F14", "HFC23", "HFC32", "HFC43_10", "HFC125", "HFC134a", "HFC143a", "HFC227ea", "HFC245fa", "SF6", "CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]
ods_names = ["CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]

set_dimension!(m, :other_ghg, other_ghg_names) # FAIR dims - Set index for Kyoto and ozone ozone-depleting gases.
set_dimension!(m, :ozone_depleting_substances, ods_names) # FAIR dims - Set index for Kyoto and ozone ozone-depleting gases.

# expand the time dimension of the model to match the full extent of FAIR and 
# explicitly the `first` attribute for the MimiFUND components
set_dimension!(m, :time, collect(FAIR_first:FAIR_last))
for comp in FUND_comp_names
    Mimi.set_firstlast!(m, comp, first = FUND_first)
end

# add the MimiFAIR components to the MimiFUND model after the emissions 
# component, noting they will take the first and last of the model which is now
# FUND_first and FUND_last
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
add_comp!(m, aerosol_indirect_rf; after = :aerosol_direct_rf)
add_comp!(m, bc_snow_rf; after = :aerosol_indirect_rf);
add_comp!(m, landuse_rf; after = :bc_snow_rf);
add_comp!(m, contrails_rf; after = :landuse_rf);
add_comp!(m, total_rf; after = :contrails_rf);
add_comp!(m, temperature; after = :total_rf);

# update the FUND component parameters with padded versions
update_MimiFUND_params!(m)

# set all the FAIR component parameters and make their internal connections


# connect FUND and FAIR

# destination: FAIR component :co2_cycle --> parameter :E_CO₂
# previous source: exogenous parameter
# new source: FUND :emissions --> variable :mco2 
rcp_emissions, volcano_forcing, solar_forcing, gas_data, gas_fractions, conversions = MimiFAIR.load_fair_data(FAIR_first, FAIR_last, rcp);
FAIR_CO₂_backup = (rcp_emissions.FossilCO2 .+ rcp_emissions.OtherCO2) * 1000
connect_param!(m, :co2_cycle, :E_CO₂, :emissions, :mco2, FAIR_CO₂_backup, backup_offset = 1)

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
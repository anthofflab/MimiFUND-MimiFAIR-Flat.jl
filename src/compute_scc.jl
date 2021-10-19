##
## COMPUTE FUNDFAIR16 SCGHGS
##


# include(joinpath(@__DIR__, "helper.jl")) # path relative to this file
include(joinpath(@__DIR__, "src/helper.jl"))

# include(joinpath(@__DIR__, "marginaldamages_MimiFAIR162.jl")) # path relative to this file
include(joinpath(@__DIR__, "src/marginaldamages_MimiFAIR162.jl"))

# include(joinpath(@__DIR__, "MimiFUND-MimiFAIR162-Flat.jl")) # path relative to this file
include(joinpath(@__DIR__, "src/MimiFUND-MimiFAIR162-Flat.jl"))


scenarios = ["ssp245", "ssp370"]
years =[2020, 2030, 2040]
rates = [0.02, 0.03]
gases = [:CO2, :CH4, :N2O]

scghg_dict = Dict()
    
for scen in scenarios
    m = get_fundfair16(ar6_scenario = scen)
    for gas in gases
        for year in years
            for rate in rates
                mm = create_marginal_fundfair16_model(m; year = year, gas = gas)
                scghg_dict[gas, year, rate, scen] = compute_fundfair16_sc(mm, year = year, rate = rate)
            end
        end
    end
end


scghg_dict


## save
using CSV # sorry david
path = "C:/Users/TTAN/Environmental Protection Agency (EPA)/NCEE Social Cost of Carbon - General/models/Notes/Code/output/fair1.6"
CSV.write(joinpath(path,"fund39_fair16_scghgs.csv"), scghg_dict)
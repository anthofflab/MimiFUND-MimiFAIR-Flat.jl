#-------------------------------------------------------
#  set up FUND and FAIR models
#-------------------------------------------------------

using MimiFAIR
using MimiFUND
using Mimi

include(joinpath(@__DIR__, "helper.jl")) # path relative to this file

# set some constants
year = 2030
pulse_size = 1e7

plots_output_dir = joinpath(@__DIR__, "..", "output")
isdir(plots_output_dir) || mkpath(plots_output_dir)

# start with default MimiFUND model
m = get_fundfair(;rcp_scenario = "RCP85", FAIR_first = 1765, FAIR_last = 2500, FUND_first = 1950, FUND_last = 2500)
run(m)

#-------------------------------------------------------
# create marginal model and perturb
#-------------------------------------------------------

pulse_size = 1e7
mm = Mimi.create_marginal_model(m, pulse_size)

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

# Add additional emissions to mm.modified
add_comp!(mm.modified, emissionspulse, after = :emissions)
nyears = length(Mimi.time_labels(m))
addem = zeros(nyears) 
baseyear = 1765

# pulse is spread over ten years, and emissions components is in Mt so divide by 1e7, and convert from CO2 to C because emissions component is in MtC
addem[(year - baseyear + 1):(year - baseyear + 1) + 9] .= pulse_size / 1e7 * 12/44

set_param!(mm.modified, :emissionspulse, :add, addem)

# Reconnect the appropriate emissions in m
connect_param!(mm.modified, :emissionspulse, :input, :emissions, :mco2, zeros(nyears), backup_offset = 1)
connect_param!(mm.modified, :multiplier, :input, :emissionspulse, :output)
run(mm)

# check results
# mm.modified
# mm.modified[:emissions, :mco2] # 736 element array, missings for pre-1950 years
# mm.modified[:emissionspulse, :input] # 736 element array, 0s for pre-1950 years
# mm.modified[:emissionspulse, :input][266:275]
# mm.modified[:emissionspulse, :add] # 736 element array, all zeros except 10 x 0.2727 between years corresponding to 2030 and 2040 (mm.modified[:emissionspulse, :add][266:275])
# mm.modified[:emissionspulse, :add][266:275]
# mm.modified[:emissionspulse, :output] # 736 element array, 0s for pre-1950 years
# mm.modified[:emissionspulse, :output][266:275]

# mm.modified[:emissionspulse, :output][266:275] == mm.modified[:emissionspulse, :input][266:275] + mm.modified[:emissionspulse, :add][266:275] # true

# mm.modified[:multiplier, :input]
# mm.modified[:multiplier, :input] == mm.modified[:emissionspulse, :output] # true
# mm.modified[:multiplier, :multiply]
# mm.modified[:multiplier, :output] # missings for pre-1950 years?
# mm.modified[:multiplier, :output][187:end] == mm.modified[:multiplier, :input][187:end] .* 0.001 # true
# mm.modified[:multiplier, :input][187:end]

# mm.modified[:co2_cycle, :E_CO₂]
# mm.modified[:co2_cycle, :E_CO₂][187:end] == mm.modified[:multiplier, :output][187:end] # true
# mm.modified[:co2_cycle, :E_CO₂][1:186] == FAIR_CO₂_backup[1:186]

#-------------------------------------------------------
# calculate SCC
#-------------------------------------------------------

marginaldamage = mm[:impactaggregation, :loss]
# marginaldamage = mm.modified[:impactaggregation, :loss] - mm.base[:impactaggregation, :loss]

# calculate discount factors
prtp = 0.05
new_years = collect(1765:1:2500)
year_index = findfirst(isequal(year), new_years)

df = zeros(length(new_years), 16)
for i in 1:length(new_years)
    if i >= year_index
        df[i,:] .= 1/(1+prtp)^(i-year_index)
    end
end

# calculate SCC
scc = sum(skipmissing(marginaldamage .* df))


## compare with FUND SCC
mfund = MimiFUND.get_model()
run(mfund)
scc_fund = MimiFUND.compute_scco2(mfund, year = 2030, eta = 0., prtp = 0.05, equity_weights = false) #...wut


#-------------------------------------------------------
#  plot results
#-------------------------------------------------------

using VegaLite
using DataFrames

# get marginal FUND model
mmfund = MimiFUND.get_marginal_model(mfund; year = year, gas = :CO2, pulse_size = pulse_size)
run(mmfund)

# get marginal FAIR model
mfair = MimiFAIR.get_model()
run(mfair)

# plot marginal damages
fundfair_md = sum(mm[:impactaggregation, :loss], dims = 2)[:,1]; 
# fundfairvals_perturbed = sum(mm.modified[:impactaggregation, :loss], dims = 2)[:,1];
# fundfairvals_base = sum(mm.base[:impactaggregation, :loss], dims = 2)[:,1];

fund_md = vcat(
    fill(missing, length(FAIR_first:FUND_first)-1),
    sum(mmfund[:impactaggregation, :loss][1:length(FUND_first:FAIR_last),:], dims = 2)
)
fund_md = fund_md[:,1]; # make a vector

df = DataFrame(
    :Year => Mimi.time_labels(m),
    :FUND => fund_md, 
    :FUNDFAIR => fundfair_md)
    # :FUNDFAIR_base => fundfairvals_base,
    # :FUNDFAIR_perturbed => fundfairvals_perturbed)

# stack(df, [:FUND, :FUNDFAIR_base, :FUNDFAIR_perturbed]) |> 
stack(df, [:FUND, :FUNDFAIR]) |> 

@vlplot(
    :line, 
    x = :Year,
    y = :value,
    color = :variable
) |>
save("fundfair_marginaldamages.png")

# plot temperature delta
fundfair_tempdelta =  mm[:temperature, :T]
fund_tempdelta = vcat(
    fill(missing, length(FAIR_first:FUND_first)-1),
    mmfund[:climateco2cycle, :temp][1:length(FUND_first:FAIR_last),:]
)
fund_tempdelta = fund_tempdelta[:,1]; # make a vector

df = DataFrame(
    :Year => Mimi.time_labels(m),
    :FUND => fund_tempdelta,
    # :FAIR => fairvals,
    :FUNDFAIR => fundfair_tempdelta
)
stack(df, [:FUND, :FUNDFAIR]) |> 
# stack(df, [:FUND, :FAIR, :FUNDFAIR]) |> 

@vlplot(
    :line, 
    x = :Year,
    y = :value,
    color = :variable
) |>
save("fundfair_tempdelta.png")

## using plots instead of vegalite

using Plots

# temp deltas for 1GtCO2 perturbation
fund_tempdelta[186:736]
fundfair_tempdelta[186:736]

# years = Mimi.time_labels(m)
years = collect(1950:1:2500)

Plots.plot(years, fund_tempdelta[186:736], 
        title = "Temperature Deltas from 1tCO2 Emissions Pulse",
        xaxis = "Year",
        ylabel = "Degrees C",
        label = "FUND")
Plots.plot!(years, fundfair_tempdelta[186:736], label = "FUNDFAIR (RCP8.5)") 

Plots.savefig(joinpath(plots_output_dir, "fundfair_tempdelta"))

# marginal damages

Plots.plot(years, fund_md[186:736], 
        title = "Marginal Damages from 1tCO2 Emissions Pulse",
        xaxis = "Year",
        ylabel = "\$",
        label = "FUND")
Plots.plot!(years, fundfair_md[186:736], label = "FUNDFAIR (RCP8.5)")

Plots.savefig(joinpath(plots_output_dir, "fundfair_marginaldamages"))

# emissions
fundfair_base_emissions = mm.base[:co2_cycle, :E_CO₂]
fund_base_emissions = mfund[:emissions, :mco2][1:551]/1000
fair_emissions = mfair[:co2_cycle, :E_CO₂][186:736]

Plots.plot(years, fundfair_base_emissions[186:736], 
        title = "CO2 Emissions",
        xaxis = "Year",
        ylabel = "GtC / year",
        label = "FUNDFAIR / FUND")
# Plots.plot!(years, fund_base_emissions, label = "FUND Base")
Plots.plot!(years, fair_emissions, label = "FAIR (RCP8.5)")

Plots.savefig("fundfair_emissions")

# atmospheric CO2
fund_acco2 = mfund[:climateco2cycle, :acco2][1:551]
fundfair_acco2 = mm.base[:co2_cycle, :C][186:736]
fair_acco2 = mfair[:co2_cycle, :C][186:736]

Plots.plot(years, fund_acco2, 
        title = "Atmospheric CO2 Concentration",
        xaxis = "Year",
        ylabel = "ppm",
        label = "FUND")
Plots.plot!(years, fundfair_acco2, label = "FUNDFAIR")
Plots.plot!(years, fair_acco2, label = "FAIR (RCP8.5)")

Plots.savefig("fundfair_acco2")

# radiative forcing

mfair85 = MimiFAIR.get_model(rcp_scenario = "RCP85")
run(mfair85)

mfair45 = MimiFAIR.get_model(rcp_scenario = "RCP45")
run(mfair45)


mfair26 = MimiFAIR.get_model(rcp_scenario = "RCP26")
run(mfair26)

fund_rf = mfund[:climateforcing, :radforc][1:551]
fundfair_rf_base = mm.base[:temperature, :F][186:736]
fundfair_rf_perturbed = mm.modified[:temperature, :F][186:736]

fair85_rf = mfair85[:temperature, :F][186:736]
fair45_rf = mfair45[:temperature, :F][186:736]
fair26_rf = mfair26[:temperature, :F][186:736]

Plots.plot(years, fair85_rf, 
        title = "Total Radiative Forcing",
        xaxis = "Year",
        ylabel = "",
        label = "FAIR (RCP8.5)",
        legend = :topleft)
Plots.plot!(years, fair45_rf, label = "FAIR (RCP4.5)")
Plots.plot!(years, fair26_rf, label = "FAIR (RCP2.6)")
Plots.plot!(years, fundfair_rf_base, label = "FUNDFAIR (Base)")
Plots.plot!(years, fundfair_rf_perturbed, label = "FUNDFAIR (Perturbed)")
Plots.plot!(years, fund_rf, label = "FUND")


Plots.savefig(joinpath(plots_output_dir, "fundfair_acco2"))


# change in rf
fundfair_rf_delta = mm[:temperature, :F][186:736]
fund_rf_delta = mmfund[:climateforcing, :radforc][1:551]

Plots.plot(years, fund_rf_delta, 
        title = "Change in Radiative Forcing from 1tCO2 Pulse",
        xaxis = "Year",
        ylabel = "",
        label = "FUND")
        # legend = :topright)
Plots.plot!(years, fundfair_rf_delta, label = "FUNDFAIR")

Plots.savefig(joinpath(plots_output_dir, "fundfair_rfdelta"))

# change in atmospheric CO2
fund_acco2_delta = mmfund[:climateco2cycle, :acco2][1:551]
fundfair_acco2_delta = mm[:co2_cycle, :C][186:736]

Plots.plot(years, fund_acco2_delta, 
        title = "Change in Atmospheric CO2 Concentration",
        xaxis = "Year",
        ylabel = "ppm",
        label = "FUND")
Plots.plot!(years, fundfair_acco2_delta, label = "FUNDFAIR")

Plots.savefig(joinpath(plots_output_dir, "fundfair_acco2_delta"))

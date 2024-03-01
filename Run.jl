
# the Revise package allows you to rerun code after making variable/function/method changes without re-compiling everything
using Revise, Tables,CSV
using Debugger
using RCall

@rlibrary seacarb

break_on(:error)

_revise_mode__ = :eval

include("AlgaeRiverReactor.jl")
using .AlgaeRiverReactor


productivity_v = zeros(8,3)
vol_fl_v = [0.800,0.700,0.600,0.500,0.400,0.300,0.200,0.100]
len_v = [15, 30, 45, 60, 75, 90, 105, 120] 
co2_v = [250, 500, 750, 1000, 1250, 1500, 1750, 2000]

for l = 1:length(vol_fl_v)
    Revise.track("AlgaeRiverReactor.jl")
    Revise.track("LoadParameters.jl")
    Revise.track("PDE_AlgaeBiomass.jl")
    Revise.track("PDE_Temperature.jl")
    Revise.track("PDE_CO2.jl")
    Revise.track("MakePlots.jl")

    @show Revise.watched_files

    productivity_v[l,1] = AlgaeRiverReactor.Run(l, 1, 1, 1)
    @show argmax(productivity_v[:,1])
end

idx = argmax(productivity_v[:,1])

for q = 1:length(len_v)
    Revise.track("AlgaeRiverReactor.jl")
    Revise.track("LoadParameters.jl")
    Revise.track("PDE_AlgaeBiomass.jl")
    Revise.track("PDE_Temperature.jl")
    Revise.track("PDE_CO2.jl")
    Revise.track("MakePlots.jl")

    @show Revise.watched_files

    productivity_v[q,2] = AlgaeRiverReactor.Run(idx, q, 1, 2)
end

idx1 = argmax(productivity_v[:,2])

for q = 1:length(co2_v)
    Revise.track("AlgaeRiverReactor.jl")
    Revise.track("LoadParameters.jl")
    Revise.track("PDE_AlgaeBiomass.jl")
    Revise.track("PDE_Temperature.jl")
    Revise.track("PDE_CO2.jl")
    Revise.track("MakePlots.jl")

    @show Revise.watched_files

    productivity_v[q,3] = AlgaeRiverReactor.Run(idx, idx1, q, 3)
    @show productivity_v
end

 
vol_v_prod_table = hcat(productivity_v[1:length(vol_fl_v),1], vol_fl_v)
CSV.write("vol_v_prod.csv", Tables.table(vol_v_prod_table), writeheader = false)

len_v_prod_table = hcat(productivity_v[1:length(len_v),2], len_v)
CSV.write("len_v_prod.csv", Tables.table(len_v_prod_table), writeheader = false)

co2_v_prod_table = hcat(productivity_v[1:length(co2_v),3], co2_v)
CSV.write("co2_v_prod.csv", Tables.table(co2_v_prod_table), writeheader = false)






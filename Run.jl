
# the Revise package allows you to rerun code after making variable/function/method changes without re-compiling everything
using Revise
using Debugger

break_on(:error)

_revise_mode__ = :eval

include("AlgaeRiverReactor.jl")
using .AlgaeRiverReactor
productivity_v = zeros(8,3)
vol_fl_v = [0.005; 0.010; 0.015; 0.020; 0.025; 0.030; 0.035; 0.040]
len_v = [15; 30; 45; 60; 75; 90; 105; 120]
co2_v = [0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8]
for l = 1:8
    Revise.track("AlgaeRiverReactor.jl")
    Revise.track("LoadParameters.jl")
    Revise.track("PDE_AlgaeBiomass.jl")
    Revise.track("PDE_Temperature.jl")
    Revise.track("PDE_CO2.jl")
    Revise.track("MakePlots.jl")

    @show Revise.watched_files

    productivity_v[l,1] = AlgaeRiverReactor.Run(l, 5, 5, 1)
    @show argmax(productivity_v[:,1])
end

idx = argmax(productivity_v[:,1])

for q = 1:8
    Revise.track("AlgaeRiverReactor.jl")
    Revise.track("LoadParameters.jl")
    Revise.track("PDE_AlgaeBiomass.jl")
    Revise.track("PDE_Temperature.jl")
    Revise.track("PDE_CO2.jl")
    Revise.track("MakePlots.jl")

    @show Revise.watched_files

    productivity_v[q,2] = AlgaeRiverReactor.Run(idx, q, 5, 2)
end

idx1 = argmax(productivity_v[:,2])

for t = 1:8
    Revise.track("AlgaeRiverReactor.jl")
    Revise.track("LoadParameters.jl")
    Revise.track("PDE_AlgaeBiomass.jl")
    Revise.track("PDE_Temperature.jl")
    Revise.track("PDE_CO2.jl")
    Revise.track("MakePlots.jl")

    @show Revise.watched_files

    productivity_v[q,3] = AlgaeRiverReactor.Run(idx, idx1, t, 3)
    @show productivity_v
end

 
vol_v_prod_table = hcat(productivity_v[:,1], vol_fl_v)
CSV.write("vol_v_prod.csv", Tables.table(vol_v_prod_table), writeheader = false)

len_v_prod_table = hcat(productivity_v[:,2], len_v)
CSV.write("len_v_prod.csv", Tables.table(len_v_prod_table), writeheader = false)

co2_v_prod_table = hcat(productivity_v[:,3], co2_v)
CSV.write("co2_v_prod.csv", Tables.table(co2_v_prod_table), writeheader = false)






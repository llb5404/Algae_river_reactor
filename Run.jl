
# the Revise package allows you to rerun code after making variable/function/method changes without re-compiling everything
using Revise, Tables,CSV
using Debugger


break_on(:error)

_revise_mode__ = :eval

include("AlgaeRiverReactor.jl")
include("LoadParameters.jl") 
using .AlgaeRiverReactor     



## Load Parameters
filesuffix1 = ["v_v1", "v_v2", "v_v3", "v_v4", "v_v5", "v_v6", "v_v7", "v_v8"]

productivity_v = zeros(8,2)
w_loss = zeros(8,2)

params = LoadDefaultParameters(filesuffix1[1], 8, 1) #constant params


for l = 1:8
    Revise.track("AlgaeRiverReactor.jl")
    Revise.track("LoadParameters.jl")
    Revise.track("PDE_AlgaeBiomass.jl")
    Revise.track("PDE_Temperature.jl")
    Revise.track("PDE_CO2.jl")
    Revise.track("MakePlots.jl")

    @show Revise.watched_files
    params = LoadDefaultParameters(filesuffix1[l], 8, 1)

    (productivity_v[l,1], w_loss[l,1]) = AlgaeRiverReactor.Run(l, 8, 1)
    

    @show productivity_v
    @show w_loss


end

idx = 1
#argmax(productivity_v[:,1])

for q = 1:8
    Revise.track("AlgaeRiverReactor.jl")
    Revise.track("LoadParameters.jl")
    Revise.track("PDE_AlgaeBiomass.jl")
    Revise.track("PDE_Temperature.jl")
    Revise.track("PDE_CO2.jl")
    Revise.track("MakePlots.jl")

    @show Revise.watched_files
    params = LoadDefaultParameters(filesuffix1[q], idx, q)

    (productivity_v[q,2], w_loss[q,2]) = AlgaeRiverReactor.Run(idx, q, 2)

    @show productivity_v
    @show w_loss
end


 
vol_v_prod_table = hcat(productivity_v[1:8,1], w_loss[1:8,1])
CSV.write("vol_v_prod.csv", Tables.table(vol_v_prod_table), writeheader = ["ACP","Water Loss (kg)"])

len_v_prod_table = hcat(productivity_v[1:8,2], w_loss[1:8,2])
CSV.write("len_v_prod.csv", Tables.table(len_v_prod_table), writeheader = ["ACP","Water Loss (kg)"])







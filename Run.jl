
# the Revise package allows you to rerun code after making variable/function/method changes without re-compiling everything
using Revise, Tables,CSV
using Debugger


using Logging
using Base.StackTraces: stacktrace, StackFrame

# Open the log file


_revise_mode__ = :eval

include("AlgaeRiverReactor.jl")
include("LoadParameters.jl") 
using .AlgaeRiverReactor     



## Load Parameters
filesuffix1 = ["v_v1", "v_v2", "v_v3", "v_v4", "v_v5", "v_v6", "v_v7", "v_v8"]

productivity_v = zeros(8,3)
co2_ratio_v = zeros(8,3)
DIC_v = zeros(8,3)
pH_out = zeros(8,3)

vol_fl_v = zeros(8,3)
len_v = zeros(8,3)
co2_vec = zeros(8,3)
final_height = zeros(8,3)
bm_o_vec = zeros(8,3)

params = LoadDefaultParameters(filesuffix1[1], 8, 1, 1) #constant params

co2_i_vec = params.co2_init*ones(8,1)
bm_i_vec = params.input_biomass_concentration*ones(8,1)
ll_i_vec = params.reactor_initial_liquid_level*ones(8,1)
mol_frac = params.mol_frac_co2*ones(8,1)
pH_in = params.pH_init*ones(8,1)



#initial co2 conc, initial bm density, initial liquid level



for l = 1:8
    Revise.track("AlgaeRiverReactor.jl")
    Revise.track("LoadParameters.jl")
    Revise.track("PDE_AlgaeBiomass.jl")
    Revise.track("PDE_Temperature.jl")
    Revise.track("PDE_CO2.jl")
    Revise.track("MakePlots.jl")

    @show Revise.watched_files
    params = LoadDefaultParameters(filesuffix1[l], 8, 1, 1)

    (productivity_v[l,1], co2_ratio_v[l,1], DIC_v[l,1],final_height[l,1],pH_out[l,1], bm_o_vec[l,1]) = AlgaeRiverReactor.Run(l, 8, 1, 1)
   
    vol_fl_v[l,1] = params.flow_rates[l]
    len_v[l,1] = params.lengths[8]
    co2_vec[l,1] = params.co2_v[1]

    @show productivity_v
    @show co2_ratio_v
    @show DIC_v
    @show vol_fl_v
    @show co2_vec
    @show len_v
    @show final_height
    @show pH_out

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
    params = LoadDefaultParameters(filesuffix1[q], idx, q, 1)

    (productivity_v[q,2], co2_ratio_v[q,2], DIC_v[q,2],final_height[q,2],pH_out[q,2], bm_o_vec[q,2]) = AlgaeRiverReactor.Run(idx, q, 1, 2)
    

    vol_fl_v[q,2] = params.flow_rates[idx]
    len_v[q,2] = params.lengths[q]
    co2_vec[q,2] = params.co2_v[1]

    @show productivity_v
    @show co2_ratio_v
    @show DIC_v
    @show vol_fl_v
    @show co2_vec
    @show len_v
    @show final_height
    @show pH_out
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

    params = LoadDefaultParameters(filesuffix1[t], idx, idx1, t)

    (productivity_v[t,3], co2_ratio_v[t,3], DIC_v[t,3],final_height[t,3],pH_out[t,3], bm_o_vec[t,3]) = AlgaeRiverReactor.Run(idx, idx1, t, 3)
   
    vol_fl_v[t,3] = params.flow_rates[idx]
    len_v[t,3] = params.lengths[idx1]
    co2_vec[t,3] = params.co2_v[t]

    @show productivity_v
    @show co2_ratio_v
    @show DIC_v
    @show vol_fl_v
    @show co2_vec
    @show len_v
    @show final_height
    @show pH_out

end

 
vol_v_prod_table = hcat(vol_fl_v[1:8,1],len_v[1:8,1], mol_frac[1:8,1], co2_vec[1:8,1], co2_i_vec[1:8,1], bm_i_vec[1:8,1], bm_o_vec[1:8,1], ll_i_vec[1:8,1], final_height[1:8,1],productivity_v[1:8,1], co2_ratio_v[1:8,1], DIC_v[1:8,1], pH_in[1:8,1], pH_out[1:8,1])
CSV.write("vol_v_prod.csv", Tables.table(vol_v_prod_table), writeheader = ["Vol Flr (m3/hr)", "Len (m)", "Mol Frac CO2", "% Flr OC", "CO2 init (g/m3)", "BM init (g/L)", "BM out (g/L)","H init (m)", "H final (m)", "ACP","CO2/DIC","DIC (uM)","In pH", "Out pH"])

len_v_prod_table = hcat(vol_fl_v[1:8,2],len_v[1:8,2], mol_frac[1:8,1], co2_vec[1:8,2], co2_i_vec[1:8,1], bm_i_vec[1:8,1], bm_o_vec[1:8,2], ll_i_vec[1:8,1], final_height[1:8,2],productivity_v[1:8,2], co2_ratio_v[1:8,2], DIC_v[1:8,2], pH_in[1:8,1], pH_out[1:8,2])
CSV.write("len_v_prod.csv", Tables.table(len_v_prod_table), writeheader = ["Vol Flr (m3/hr)", "Len (m)", "Mol Frac CO2", "% Flr OC", "CO2 init (g/m3)", "BM init (g/L)", "BM out (g/L)", "H init (m)", "H final (m)", "ACP","CO2/DIC","DIC (uM)", "In pH", "Out pH"])

co2_v_prod_table = hcat(vol_fl_v[1:8,3],len_v[1:8,3], mol_frac[1:8,1], co2_vec[1:8,3], co2_i_vec[1:8,1], bm_i_vec[1:8,1], bm_o_vec[1:8,3], ll_i_vec[1:8,1], final_height[1:8,3],productivity_v[1:8,3], co2_ratio_v[1:8,3], DIC_v[1:8,3], pH_in[1:8,1], pH_out[1:8,3])
CSV.write("co2_v_prod.csv", Tables.table(co2_v_prod_table), writeheader = ["Vol Flr (m3/hr)", "Len (m)", "Mol Frac CO2", "% Flr OC", "CO2 init (g/m3)", "BM init (g/L)", "BM out (g/L)", "H init (m)", "H final (m)", "ACP","CO2/DIC","DIC (uM)", "In pH", "Out pH"])






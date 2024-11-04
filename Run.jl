
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

productivity_v = zeros(8,3)
sp_growth_v = zeros(8,3)
co2_ratio_v = zeros(8,3)
DIC_v = zeros(8,3)
pH_out = zeros(8,3)

vol_fl_v = zeros(8,3)
len_v = zeros(8,3)
strip_conc = zeros(8,3)
co2_vec = zeros(8,3)
final_height = zeros(8,3)
bm_o_vec = zeros(8,3)
w_loss = zeros(8,3)

params = LoadDefaultParameters(filesuffix1[1], 8, 1, 1) #constant params

co2_i_vec = (params.co2_init/1000)*ones(8,1)
bm_i_vec = params.input_biomass_concentration*ones(8,1)
ll_i_vec = params.reactor_initial_liquid_level*ones(8,1)
mol_frac = params.mol_frac_co2*ones(8,1)
pH_in = params.pH_init*ones(8,1)
strip_conc_vec = (params.strip_concentration/1000)*ones(8,1)



#initial co2 conc, initial bm density, initial liquid level



for l = 1:8
    Revise.track("AlgaeRiverReactor.jl")
    Revise.track("LoadParameters.jl")
    Revise.track("PDE_AlgaeBiomass.jl")
    Revise.track("PDE_Temperature.jl")
    Revise.track("PDE_CO2.jl")
    Revise.track("MakePlots.jl")

    @show Revise.watched_files
    params = LoadDefaultParameters(filesuffix1[l], 1, 1,1)

    ((productivity_v[l,1], co2_ratio_v[l,1], DIC_v[l,1],final_height[l,1],pH_out[l,1], bm_o_vec[l,1], sp_growth_v[l,1]), w_loss[l,1]) = AlgaeRiverReactor.Run(l, 1, 1,1)
    
    vol_fl_v[l,1] = params.flow_rates[l]
    len_v[l,1] = params.lengths[1]
    co2_vec[l,1] = params.strip_flr[1]
    strip_conc[l,1] = params.s_conc[1]/1000

    @show productivity_v
    @show co2_ratio_v
    @show DIC_v
    @show vol_fl_v
    @show co2_vec
    @show len_v
    @show final_height
    @show pH_out
    @show w_loss
    @show sp_growth_v

end

idx = 8

for q = 1:8
    Revise.track("AlgaeRiverReactor.jl")
    Revise.track("LoadParameters.jl")
    Revise.track("PDE_AlgaeBiomass.jl")
    Revise.track("PDE_Temperature.jl")
    Revise.track("PDE_CO2.jl")
    Revise.track("MakePlots.jl")

    @show Revise.watched_files
    params = LoadDefaultParameters(filesuffix1[q], idx, q, 1)

    ((productivity_v[q,2], co2_ratio_v[q,2], DIC_v[q,2],final_height[q,2],pH_out[q,2], bm_o_vec[q,2], sp_growth_v[q,2]), w_loss[q,2]) = AlgaeRiverReactor.Run(idx, q ,1,2)
    

    vol_fl_v[q,2] = params.flow_rates[idx]
    len_v[q,2] = params.lengths[q]
    strip_conc[q,2] = params.s_conc[q]/1000
    co2_vec[q,2] = params.strip_flr[1]

    @show productivity_v
    @show co2_ratio_v
    @show DIC_v
    @show vol_fl_v
    @show co2_vec
    @show len_v
    @show final_height
    @show pH_out
    @show w_loss
    @show sp_growth_v
end

idx1 = 1

for t = 1:8
    Revise.track("AlgaeRiverReactor.jl")
    Revise.track("LoadParameters.jl")
    Revise.track("PDE_AlgaeBiomass.jl")
    Revise.track("PDE_Temperature.jl")
    Revise.track("PDE_CO2.jl")
    Revise.track("MakePlots.jl")

    @show Revise.watched_files

    params = LoadDefaultParameters(filesuffix1[t], idx, idx1, t)

    ((productivity_v[t,3], co2_ratio_v[t,3], DIC_v[t,3],final_height[t,3],pH_out[t,3], bm_o_vec[t,3], sp_growth_v[t,3]), w_loss[t,3]) = AlgaeRiverReactor.Run(idx, idx1, t, 3)
   
    vol_fl_v[t,3] = params.flow_rates[idx]
    len_v[t,3] = params.lengths[idx1]
    strip_conc[t,3] = params.s_conc[idx1]/1000
    co2_vec[t,3] = params.strip_flr[t]

    @show productivity_v
    @show co2_ratio_v
    @show DIC_v
    @show vol_fl_v
    @show co2_vec
    @show len_v
    @show final_height
    @show pH_out
    @show w_loss
    @show sp_growth_v

end


 
vol_v_prod_table = hcat(vol_fl_v[1:8,1],len_v[1:8,1], mol_frac[1:8,1], co2_vec[1:8,1], co2_i_vec[1:8,1], bm_i_vec[1:8,1], bm_o_vec[1:8,1], ll_i_vec[1:8,1], final_height[1:8,1],productivity_v[1:8,1], co2_ratio_v[1:8,1], DIC_v[1:8,1], pH_in[1:8,1], pH_out[1:8,1], strip_conc_vec[1:8,1], w_loss[1:8,1], sp_growth_v[1:8,1], strip_conc[1:8,1])
CSV.write("vol_v_prod.csv", Tables.table(vol_v_prod_table), writeheader = ["Vol Flr (m3/hr)", "Len (m)", "Mol Frac CO2", "Strip Column Flr (gpm)", "CO2 init (g/L)", "BM init (g/L)", "BM out (g/L)","H init (m)", "H final (m)", "ACP","CO2/DIC","DIC (uM)","In pH", "Out pH", "Strip Conc (g/L)", "Water Loss (kg)", "Specific Growth Rate (hr-1)", "Strip Conc (g/L)"])

len_v_prod_table = hcat(vol_fl_v[1:8,2],len_v[1:8,2], mol_frac[1:8,1], co2_vec[1:8,2], co2_i_vec[1:8,1], bm_i_vec[1:8,1], bm_o_vec[1:8,2], ll_i_vec[1:8,1], final_height[1:8,2],productivity_v[1:8,2], co2_ratio_v[1:8,2], DIC_v[1:8,2], pH_in[1:8,1], pH_out[1:8,2], strip_conc_vec[1:8,1], w_loss[1:8,2], sp_growth_v[1:8,2], strip_conc[1:8,2])
CSV.write("len_v_prod.csv", Tables.table(len_v_prod_table), writeheader = ["Vol Flr (m3/hr)", "Len (m)", "Mol Frac CO2", "Strip Column Flr (gpm)", "CO2 init (g/L)", "BM init (g/L)", "BM out (g/L)", "H init (m)", "H final (m)", "ACP","CO2/DIC","DIC (uM)", "In pH", "Out pH", "Strip Conc (g/L)", "Water Loss (kg)", "Specific Growth Rate (hr-1)", "Strip Conc (g/L)"])

co2_v_prod_table = hcat(vol_fl_v[1:8,3],len_v[1:8,3], mol_frac[1:8,1], co2_vec[1:8,3], co2_i_vec[1:8,1], bm_i_vec[1:8,1], bm_o_vec[1:8,3], ll_i_vec[1:8,1], final_height[1:8,3],productivity_v[1:8,3], co2_ratio_v[1:8,3], DIC_v[1:8,3], pH_in[1:8,1], pH_out[1:8,3], strip_conc_vec[1:8,1], w_loss[1:8,3], sp_growth_v[1:8,3], strip_conc[1:8,3])
CSV.write("co2_v_prod.csv", Tables.table(co2_v_prod_table), writeheader = ["Vol Flr (m3/hr)", "Len (m)", "Mol Frac CO2", "Strip Column Flr (gpm)", "CO2 init (g/L)", "BM init (g/L)", "BM out (g/L)", "H init (m)", "H final (m)", "ACP","CO2/DIC","DIC (uM)", "In pH", "Out pH", "Strip Conc (g/L)", "Water Loss (kg)", "Specific Growth Rate (hr-1)", "Strip Conc (g/L)"])





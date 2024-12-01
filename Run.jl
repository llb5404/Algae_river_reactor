
# the Revise package allows you to rerun code after making variable/function/method changes without re-compiling everything
using Revise, Tables,CSV
using Debugger

_revise_mode__ = :eval

include("AlgaeRiverReactor.jl")
include("LoadParameters.jl") 
using .AlgaeRiverReactor     

## Load Parameters
params = LoadDefaultParameters(8, 1) #constant params


##Initialize Vectors for Model Outputs
productivity_v = zeros(8,3)
sp_growth_v = zeros(8,3)
DIC_v = zeros(8,3)
pH_out = zeros(8,3)
final_height = zeros(8,3)
bm_o_vec = zeros(8,3)
w_loss = zeros(8,3)

#ACP, specific growth at outlet, dissolved inorganic carbon, pH at outlet, height at outlet, biomass conc at outlet, water loss across reactor

## Initialize Vectors for Model Variables
vol_fl_v = zeros(8,3)
len_v = zeros(8,3)
strip_conc = zeros(8,3)
co2_vec = zeros(8,3)

#volumetric flowrate (inlet), length, absorption column outlet concentration, absorption column flowrate

#Initialize Vectors for Model Constants
co2_i_vec = (params.co2_init/1000)*ones(8,1)
bm_i_vec = params.input_biomass_concentration*ones(8,1)
ll_i_vec = params.reactor_initial_liquid_level*ones(8,1)
pH_in = params.pH_init*ones(8,1)

#inlet co2 concentration, inlet biomass concentration, initial liquid level, initial pH

#Volumetric Flowrate

idx_q = 5 # fix the second variable (absorption column concentration)

for l = 1:8 
    Revise.track("AlgaeRiverReactor.jl")
    Revise.track("LoadParameters.jl")
    Revise.track("PDE_AlgaeBiomass.jl")
    Revise.track("PDE_Temperature.jl")
    Revise.track("PDE_CO2.jl")
    Revise.track("MakePlots.jl")

    @show Revise.watched_files
    params = LoadDefaultParameters(1, idx_q)

    (productivity_v[l,1], DIC_v[l,1],final_height[l,1],pH_out[l,1], bm_o_vec[l,1], sp_growth_v[l,1], w_loss[l,1]) = AlgaeRiverReactor.Run(l, idx_q)
    
    vol_fl_v[l,1] = params.flow_rates[l]
    len_v[l,1] = params.lengths[1]
    co2_vec[l,1] = params.volumetric_flow_rate_strip
    strip_conc[l,1] = params.s_conc[idx_q]/1000

    @show productivity_v
    @show DIC_v
    @show vol_fl_v
    @show co2_vec
    @show len_v
    @show final_height
    @show pH_out
    @show w_loss
    @show sp_growth_v

end

#Absorption Column Concentration 

idx_l = 1 #fix the first variable (volumetric flowrate)

for q = 1:8
    Revise.track("AlgaeRiverReactor.jl")
    Revise.track("LoadParameters.jl")
    Revise.track("PDE_AlgaeBiomass.jl")
    Revise.track("PDE_Temperature.jl")
    Revise.track("PDE_CO2.jl")
    Revise.track("MakePlots.jl")

    @show Revise.watched_files
    params = LoadDefaultParameters(idx_l, q)

    (productivity_v[q,2], DIC_v[q,2],final_height[q,2],pH_out[q,2], bm_o_vec[q,2], sp_growth_v[q,2], w_loss[q,2]) = AlgaeRiverReactor.Run(idx_l, q)

    vol_fl_v[q,2] = params.flow_rates[idx_l]
    len_v[q,2] = params.lengths[q]
    strip_conc[q,2] = params.s_conc[q]/1000
    co2_vec[q,2] = params.volumetric_flow_rate_strip

    @show productivity_v
    @show DIC_v
    @show vol_fl_v
    @show co2_vec
    @show len_v
    @show final_height
    @show pH_out
    @show w_loss
    @show sp_growth_v
end


#Generate Tables
vol_v_prod_table = hcat(vol_fl_v[1:8,1],len_v[1:8,1], co2_vec[1:8,1], co2_i_vec[1:8,1], bm_i_vec[1:8,1], bm_o_vec[1:8,1], ll_i_vec[1:8,1], final_height[1:8,1],productivity_v[1:8,1], DIC_v[1:8,1], pH_in[1:8,1], pH_out[1:8,1], w_loss[1:8,1], sp_growth_v[1:8,1], strip_conc[1:8,1])
CSV.write("vol_v_prod.csv", Tables.table(vol_v_prod_table), writeheader = ["Vol Flr (m3/hr)", "Len (m)", "Strip Column Flr (m3/hr)", "CO2 init (g/L)", "BM init (g/L)", "BM out (g/L)","H init (m)", "H final (m)", "ACP","DIC (uM)","In pH", "Out pH", "Water Loss (kg)", "Specific Growth Rate (hr-1)", "Strip Conc (g/L)"])

len_v_prod_table = hcat(vol_fl_v[1:8,2],len_v[1:8,2], co2_vec[1:8,2], co2_i_vec[1:8,1], bm_i_vec[1:8,1], bm_o_vec[1:8,2], ll_i_vec[1:8,1], final_height[1:8,2],productivity_v[1:8,2], DIC_v[1:8,2], pH_in[1:8,1], pH_out[1:8,2], w_loss[1:8,2], sp_growth_v[1:8,2], strip_conc[1:8,2])
CSV.write("len_v_prod.csv", Tables.table(len_v_prod_table), writeheader = ["Vol Flr (m3/hr)", "Len (m)", "Strip Column Flr (m3/hr)", "CO2 init (g/L)", "BM init (g/L)", "BM out (g/L)", "H init (m)", "H final (m)", "ACP","DIC (uM)", "In pH", "Out pH", "Water Loss (kg)", "Specific Growth Rate (hr-1)", "Strip Conc (g/L)"])






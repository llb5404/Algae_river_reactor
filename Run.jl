
# the Revise package allows you to rerun code after making variable/function/method changes without re-compiling everything
using Revise, Tables,CSV
using Debugger

_revise_mode__ = :eval

include("AlgaeRiverReactor.jl")
include("LoadParameters.jl") 
using .AlgaeRiverReactor     

## Load Parameters
params = LoadDefaultParameters(1, 1) #constant params


##Initialize Vectors for Model Outputs
productivity_v = zeros(1,1)
avg_rt = zeros(1,1)
sp_growth_v = zeros(1,1)
DIC_v = zeros(1,1)
pH_out = zeros(1,1)
final_height = zeros(1,1)
bm_o_vec = zeros(1,1)
w_loss = zeros(1,1)
avg_co2_up = zeros(1,1)

#ACP, specific growth at outlet, dissolved inorganic carbon, pH at outlet, height at outlet, biomass conc at outlet, water loss across reactor

## Initialize Vectors for Model Variables
vol_fl_v = zeros(1,1)
len_v = zeros(1,1)
strip_conc = zeros(1,1)
co2_vec = zeros(1,1)

#volumetric flowrate (inlet), length, absorption column outlet concentration, absorption column flowrate

#Initialize Vectors for Model Constants
co2_i_vec = (params.co2_init/1000)*ones(1,1)
bm_i_vec = params.input_biomass_concentration*ones(1,1)
ll_i_vec = params.reactor_initial_liquid_level*ones(1,1)
pH_in = params.pH_init*ones(1,1)

#inlet co2 concentration, inlet biomass concentration, initial liquid level, initial pH

#Volumetric Flowrate

idx_q = 1 # fix the second variable (absorption column concentration)

for l = 1:1 
    Revise.track("AlgaeRiverReactor.jl")
    Revise.track("LoadParameters.jl")
    Revise.track("PDE_AlgaeBiomass.jl")
    Revise.track("PDE_Temperature.jl")
    Revise.track("PDE_CO2.jl")
    Revise.track("MakePlots.jl")

    @show Revise.watched_files
    params = LoadDefaultParameters(1, idx_q)

    (productivity_v[l,1], DIC_v[l,1],final_height[l,1], bm_o_vec[l,1], sp_growth_v[l,1], w_loss[l,1], avg_rt[l,1], avg_co2_up[l,1]) = AlgaeRiverReactor.Run(l, idx_q)
    
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
    @show avg_rt
    @show avg_co2_up

end
#Generate Tables
vol_v_prod_table = hcat(vol_fl_v[1,1],len_v[1,1], co2_vec[1,1], co2_i_vec[1,1], bm_i_vec[1,1], bm_o_vec[1,1], ll_i_vec[1,1], final_height[1,1],productivity_v[1,1], DIC_v[1,1], pH_in[1,1], w_loss[1,1], sp_growth_v[1,1], strip_conc[1,1], avg_rt[1,1], avg_co2_up[1,1])
CSV.write("output_table.csv", Tables.table(vol_v_prod_table), writeheader = ["Vol Flr (m3/hr)", "Len (m)", "Strip Column Flr (m3/hr)", "CO2 init (g/L)", "BM init (g/L)", "BM out (g/L)","H init (m)", "H final (m)", "ACP","DIC (uM)","In pH", "Out pH", "Water Loss (kg)", "Specific Growth Rate Outlet (hr-1)", "Strip Conc (g/L)", "Average RT (hr)", "Average CO2 Uptake (g/m2/day)"])





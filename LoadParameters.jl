
## LOAD PARAMETERS
## Options:
#   Load default parameters
#   Load parameters from an Excel file

using Debugger
using Polynomials

break_on(:error)

function LoadDefaultParameters(filesuffix, l, q, t)
    ## PDE Discretization
    num_odes_y = 100
    num_odes_z = 10
    time_end = 1020       #hours
    time_interval = 1.0     #hours
    ## Physical Constants
    reference_temperature = 20.0                                        # deg Celsius
    reference_pressure = 101315.0                                       # Pascals
    gas_constant = 8.31446261815324                                     # J / K / mole
    stefan_boltzmann_constant = 5.67037442E-08                          # W / m^2 / K^4
    acceleration_gravity = 9.81*(3600^2)                                        # m/hr2
    # Air Properties
    density_air(T) = reference_pressure ./ (287.058 .* T)                 # kg/m3, T in K
    dynamic_viscosity_air(T) = 1.458E-06 .* (T.^1.5) ./ (T .+ 110.4)        # Pa sec [kg/m/sec]
    thermal_conductivity_air(T) = 0.02624 .* (T/300).^0.8646              # W/m*K
    specific_heat_capacity_air(T) = 1002.5 .+ (275E-06) .* (T .- 200).^2      # J/kg*K
    prandtl_air(T) = specific_heat_capacity_air(T) .* dynamic_viscosity_air(T) ./ thermal_conductivity_air(T)             # dimensionless
    emissivity_air = 0.80                                               # dimensionless
    molecular_weight_air = 0.02896                                      # kg/mole [average of air mixture]

    # Water Properties
    density_water_20degC = 998.2071                                             # kg / m^3
    expansion_coefficient_water = 0.0002                                        # 1/degC
    density_water(T) = density_water_20degC ./ (1 .+ expansion_coefficient_water .* (T .- (reference_temperature+273)))   #kg/m^3
    dynamic_viscosity_water(T) = (2.414E−05 .* 10^( 247.8 ./ (T .-140)))*3600         # [kg/m/hr]
    
    kinematic_viscosity_water(T) = dynamic_viscosity_water(T)/density_water(T) #m2/hr
        #Change in water properties with salinity
    A(S) = 5.328 - 9.76 * 10 ^ (-2) * S + 4.04*10^(-4)*(S)^ 2 #unitless
    B(S) = -6.913 * 10 ^ (-3) + 7.351 * 10 ^ (-4) * (S) - 3.15*10^(-6)*(S)^2 #unitless
    C(S) = 9.6 * 10 ^ (-6) - 1.927 * 10 ^ (-6) * (S) + 8.23 * 10^(-9) *(S)^2 #unitless
    D(S) = 2.5 * 10 ^ (-9) + 1.666 * 10 ^ (-9) * (S) - 7.125 * 10^(-12)*(S)^2 #unitless
    specific_heat_capacity_solution(T,S) = 1000*(A(S) + B(S)*T + C(S)*(T^2) + D(S)*(T^3)) #J/kg/K
    #above information obtained from: 
    
    #K.G. Nayar, M.H. Sharqawy, L.D. Banchik, and J.H. Lienhard V, "Thermophysical properties of seawater: A review and new correlations that include pressure dependence," Desalination, Vol. 390, pp.1-24, 2016. doi:10.1016/j.desal.2016.02.024 (preprint)
    #Mostafa H. Sharqawy, John H. Lienhard V, and Syed M. Zubair, "Thermophysical properties of seawater: A review of existing correlations and data," Desalination and Water Treatment, Vol. 16, pp.354-380, April 2010. (PDF file which includes corrections through June 2017.)

    specific_heat_capacity_water(T) = specific_heat_capacity_solution(T,0) #J/kg/K             
    thermal_conductivity_water(T) = 0.565464 .+ 1.756E-03 .* T .- 6.46E-06 .* T.^2   # J/m/K/sec
    solar_reflectance_water = 0.1357 # fraction of light reflected
    diffusion_coeff_water_air(T) = 22.5E-06 .* ( (T .+ 273.15) ./ 273.15).^1.8        # m^2 / sec
    saturated_vapor_pressure_water_air(T) = (exp(77.3450+0.0057*T-7235/T))/(T^8.2)  #Pascals
    heat_vaporization_water(T) = 2430159                                        # J / kg, T = 30 degC            [40660 * ( (647.3 - (T + 273.15)) / (647.3 - 373.2))^0.38 #J/mole, but isn't very accurate]
    emissivity_water = 0.97                                                     # dimensionless
    molecular_weight_water = 0.01801528                                         # kg/mole
    density_water_vapor(T) = 0.804 #kg/m3
    
    # CO2 Properties
    diffusion_coeff_co2_water(T) = (0.0816 .* T .- 22.601) .* 1.0E-09*3600              # m^2 / sec [equation fitted based on data in https://pubs.acs.org/doi/10.1021/je401008s]
    @show diffusion_coeff_co2_water(298)
    solubility_co2_water(T) = (0.00025989 .* (T .- 273.15).^2 .- 0.03372247 .* (T .- 273.15) .+ 1.31249383) ./ 1000.0  #mole fraction of dissolved CO2 in water at equilibrium [equation fitted based on data from https://srd.nist.gov/JPCRD/jpcrd427.pdf]
    molecular_weight_co2 = 0.04401                                              # kg / mole

    # Construction Material Properties
    thermal_diffusivity_concrete = 1011.83E0-9                                  #
    thermal_conductivity_concrete = 1.3                                         # W / m / K

    ## Geographic Properties -- Seasonal and daily variations
    global_horizontal_irradiance_data = zeros(8760,1)   # W/m^2
    ambient_temperature_data = zeros(8760,1)            # deg Celsius
    relative_humidity_data = zeros(8760,1)              # percent humidity
    wind_speed_data = zeros(8760,1)                     # m/s
    #days_data = zeros(8760,1)                           # days
    #times_data = zeros(8760,1)                          # hours

    # Reading Measured GHI Values from Pheonix, AZ [stored in "722789TYA.csv"]
    csv_reader = CSV.File("722789TYA.csv", skipto=3, header=["col1", "col2", "col3", "col4", "col5", "col6","col7", "col8", "col9", "col10", "col11", "col12","col13", "col14", "col15", "col16", "col17", "col18","col19", "col20", "col21", "col22", "col23", "col24","col25", "col26", "col27", "col28", "col29", "col30","col31", "col32", "col33", "col34", "col35", "col36","col37", "col38", "col39", "col40", "col41", "col42","col43", "col44", "col45", "col46", "col47"])
    for (i,row) in enumerate(csv_reader)
        #days_data[i] = row.col1                         # days
        #times_data[i] = row.col2                        # hours
        global_horizontal_irradiance_data[i] = convert(Float64, row.col5)     # W/m^2
        ambient_temperature_data[i] = convert(Float64, row.col32) + 273.15    # Kelvin
        relative_humidity_data[i] = convert(Float64, row.col38)               # percent humidity
        wind_speed_data[i]= convert(Float64, row.col47)                       # m/s
    end

    #April 1st: 2161
    #Nov 1st: 7297
    data_begin = 2161

    ## River Reactor Geometric Properties

    lengths = [15, 30, 45, 60, 75, 90, 105, 120]                           # meters
    reactor_length = lengths[q]

    reactor_width =  0.5                               # meters
    reactor_initial_liquid_level = 0.10             # meters
    reactor_depth = 0   # meters
    reactor_incline = (asin(reactor_depth/reactor_length))  # radians
    
    flow_rates = [0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040]
    volumetric_flow_rate_o = flow_rates[l] #m^3/hr
    

    tot_energy(H) = volumetric_flow_rate_o*(dynamic_viscosity_water(292.5))/((2*(H^2)/3)*H*reactor_width) #kg/hr^2-m2
    tot_energy_mod(H) = (volumetric_flow_rate_o)*(dynamic_viscosity_water(292.5))/((2*(H^2)/3)*H*reactor_width) #kg/hr^2-m2
    proj_react(H) = sqrt(reactor_length^2 - (tot_energy_mod(H)/(density_water_20degC*acceleration_gravity))^2) #m
    ## River Reactor Operating Parameters
    roughness_height = 5E-03 #m, taken from value for ordinary wood (need source on clay) https://www.engineeringtoolbox.com/surface-roughness-ventilation-ducts-d_209.html 
    hydraulic_diameter(H) = 4 * (reactor_width * H) / (reactor_width + 2 * H) #m
    global_shear_velocity(H) = ((hydraulic_diameter(H)*(tot_energy(H)))/density_water_20degC)^(1/2) #m/hr
    local_shear_velocity(H) = ((H*(tot_energy(H)))/density_water_20degC)^(1/2) #m/hr
    rough_reynolds_number(H,T) = (global_shear_velocity(H)*roughness_height)/(kinematic_viscosity_water(T)) #unitless
    #boundary layer heights listed in ascending order for higher Re
    boundary_layer_height_1(H,T) = roughness_height*(0.11/rough_reynolds_number(H,T)) #m
    boundary_layer_height_2(H,T) = roughness_height*(max((0.0275 −0.007*(sin((rough_reynolds_number(H,T)−4)/14)*pi)),0)^(1/2)) #m
    boundary_layer_height_3(H,T) = roughness_height*(max((0.0205 +(0.0125/sqrt(2))*(1+sin((rough_reynolds_number(H,T)−40.5/59))*pi)), 0)^(1/2)) #m
    boundary_layer_height_4 = roughness_height*0.033 #m
    
    wake_function(z,H) = (2*0.10/0.4)sin((pi/2)*((H-z)/H))^2 #unitless, used for flow at the top 80% of the reactor
    #for Turbulent Flow
        #use for j = 4Nz/5:Nz
        velocity_profile_nw(z,H,B) = global_shear_velocity(H)*((1/0.4)*log(max(((H-z)+B)/roughness_height, 1)) -2.5*log(B/roughness_height))  #m/hr
        #use for j = 0:4Nz/5
        velocity_profile_w(z,H,B) = global_shear_velocity(H)*((1/0.4)*log(max(((H-z)+B)/roughness_height, 1)) -2.5*log(B/roughness_height)+wake_function(z,H)) #m/hr
        ## Han et al., (2017) https://www.hindawi.com/journals/mpe/2018/6491501/
        ## all other info above in RR Operating Parameters from Bonakdari et al., (2008) https://www.researchgate.net/publication/225456790_Turbulent_velocity_profile_in_fully-developed_open_channel_flows  
        average_flow_velocity(H,B) = (global_shear_velocity(H)/H)*(2.5*((H+B)*log((H+B)/roughness_height)-H)-2.5*((H/5+B)*log((H/5+B)/roughness_height)-H/5)-2.5*log(B/roughness_height)*(4*H/5)+(5*(sin(pi/5)+4*pi)*H)/(20*pi)+2.5*((H/5+B)*log((H/5+B)/roughness_height)-H/5)-2.5((B)*log(B/roughness_height))-2.5*log(B/roughness_height)*(H/5)) #m/hr
    #for Laminar Flow
        velocity_profile_lam(z,H,T) = ((tot_energy(H))/(dynamic_viscosity_water(T)))*((H-z)*H +0.5*(H-z)^2) #m/hr
        average_flow_velocity_lam(H,T) = ((tot_energy(H))/(dynamic_viscosity_water(T)))*(2*(H^2)/3) #m/hr
       
    @show velocity_profile_nw(0.25,0.25,boundary_layer_height_3(0.25,298))
    @show velocity_profile_w(0,0.25,boundary_layer_height_3(0.25,298))
    volumetric_flow_rate(H,V) = reactor_width*H*V #V is Vavg
    #Reynolds number for a duct (open surface), Re = Vavg * density * hydraulic diameter / viscosity
    reynolds_number(H,T,V) = V .* density_water(T) .* hydraulic_diameter(H) ./ (dynamic_viscosity_water(T))
 
    
    input_hydraulic_diameter = 4 * (reactor_width * reactor_initial_liquid_level) / (reactor_width + 2 * reactor_initial_liquid_level) # meters
    input_temperature = 293.15                          # Kelvin
    ground_temperature = 290.0                          # Kelvin

    ## Biomass Properties
    input_biomass_concentration = 10.0 / 1000.0         # kg/m^3
    max_biomass_concentration = 10000.0 / 1000.0        # kg/m^3 where light is 100% absorbed
    max_biomass_light_saturation = 900.0 / 4.57         # 900 μmol m−2 s−1 converted to W/m^2
    max_biomass_specific_growth_rate = 6.75/24   # 1 / hour (obtained from Krishnan at T_opt, 0 salinity)
    
    #biomass growth adjustment factors
    salinity_factor(S) = -6E-05*S^2 + 0.0007*S + 1.0078 #unitless, S in kg/m3
    ## obtained from 2-degree polynomial fit of data obtained from Krishnan et al.
    T_opt = 308 #in C, optimum temperature
    T_min = 277.2 #in C, minimum temperature
    T_max = 320 #in C, maximum temperature
    temperature_factor_g(T) = (T_opt-T_min)*(T-T_opt) #unitless
    temperature_factor_f(T) = (T_opt-T_max)*(T_opt+T_min-2*(T)) #unitless
    temperature_factor(T) = ((T-T_max)*(T-T_min)^2)/((T_opt-T_min)*(temperature_factor_g(T)-temperature_factor_f(T))) #unitless
    pH_opt = 7.5
    pH_max = 10
    pH_min = 5
    pH_factor_g(H) = (pH_opt-pH_min)*(H-pH_opt)
    pH_factor_f(H) = (pH_opt-pH_max)*(pH_opt+pH_min-2*(H))
    pH_factor(H) = 1
    #((H-pH_max)*(H-pH_min)^2)/((pH_opt-pH_min)*(pH_factor_g(H)-pH_factor_f(H)))
    @show pH_factor(10^-14)
    


    ## opt, min, max temperatures given by Krishnan et al 
    ##temperature factor obtained from model given by Greene et al https://www.osti.gov/servlets/purl/1806213

    photosynthetic_efficiency = 0.025                   # fraction of sunlight converted to chemical energy during photosynthesis
    threshold_dissolved_co2_growth = 0.5                # kg CO2 / m^3 water, minimum dissolved CO2 concentration before growth begins to slow
    Vol(dz) = (reactor_length/num_odes_y)*dz*reactor_width #m3
    Iave(M_biomass_z, GHI, z, dz) = GHI .* 0.45 .* (1 - sum(M_biomass_z[1:z]*dz)./ max_biomass_concentration) #unitless
    phiL(M_biomass_z, GHI, z, dz) = Iave(M_biomass_z, GHI, z, dz) .* exp(1 - Iave(M_biomass_z, GHI, z, dz)./max_biomass_light_saturation) ./ max_biomass_light_saturation #unitless
    
    co2_availability_factor(C_co2) = (C_co2/1000) ./ ((C_co2/1000) + 0.00036*molecular_weight_co2 + (C_co2/1000)^2/(10*molecular_weight_co2)) #unitless, input is CO2 in g/m3
    biomass_specific_growth_rate(T, S, M_biomass_z, GHI, z, dz, C_co2,H) = max_biomass_specific_growth_rate .* phiL(M_biomass_z, GHI, z, dz) .* co2_availability_factor(C_co2).*temperature_factor(T).*salinity_factor(S)*pH_factor(H) #1/hr
    
    
    co2_per_biomass = 0.70   #based on 2.661 kg CO2 emitted when burning 1 gallon of algae (about 3.79 kg)

 
    ammonia_N_content = 0.82 #ammonia N content kg N/kg NH4
    DAP_N_content = 0.18 #di-ammonium phosphate N content kg N/kg DAP
    DAP_P_content = 0.20 #di-ammonium phosphate P content kg P/kg DAP
    N_per_biomass = 0.0450 # kg/kg algae
    P_per_biomass = 0.0120 #kg/kg algae
    #determined from theshold for algae blooms - more experimental work is needed: 

    #https://pdf.sciencedirectassets.com/271768/1-s2.0-S0043135400X0348X/1-s2.0-0043135472901820/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEBUaCXVzLWVhc3QtMSJHMEUCIGzqafx1pmG8YqjF876REGo3PwKrrVpspH2WfYa6ccuCAiEA5izjI3zuJVLWi1wAra66BdVHSCVwiqfLJuiH%2FgnHGB8qsgUILhAFGgwwNTkwMDM1NDY4NjUiDI7Qm2MoFAIfvOCSZSqPBUPoreJ7QpwQ05steBaLMuz5CaUR%2F1dHW2yHWdlCWQpSFe%2Fb5KQKUZV8cgDoCgsYUH5TMaWzOxKS6%2BveI6wt9bE9lmjzj%2BcYjHb8jK6y72mOqYDIf8Vuu41RPCPelEwxdi1RW8rSigRj6s6qzo1ksbcPVraU7USDYsNp00pzxAbMnWw8aRZs8GdX1oU9DCLfQ3h92yQBLlxJmDoam1xDhJ6lhMJ8vePRx3WJuKTLcnUR35l5N0IMcH5xKaqbvggLuCaoIKfHsfyKHh8IFTd2qh3guyIfUDViEl8mv0iERyAWcBKYkFky6VE9%2BxlvyqWR8DLuUmZHiiiWqMPaLAvvRDIlh31u4R3TxUDlDd0h17uG4ylzof1%2FFl2IfK9T%2F%2FgsIM6aaWoXz0iS38xoWvnSpk1FaISDgLpomkJo0vMXcMgppVoZFRvLm0SSMIFgmJHsxcKugU9gqze5JKnjTtu7aQ0kmehbMCfe1yIzVlc3GKsDz3LKWpZuslphmxERT%2Ff8okKXRK3JbITP%2BSINXYY6GJgqBik26fVInML1Ou0R1GbHV08MznveEozf4epcz9Wska%2BuBbSOkT5rg%2FHvrOHUInW7etCGAOZFecBvNvtGo3pNKe5LZQjxL83gcNa6aKI2G2NrHGj1aKVHUNssltmZWJgVxzGJOI4Gf1FpWAGz9Bcs7Mfe6%2Bm4Nt3ceyo8vjWuzQ%2BcokcOVqqzCLJJF39ABL7kk1F4Bii%2FgM2Hjfh1R4Osb7pIYAVpEbk9AJ0HFOmYMtchNCz71jh5koEBTzDxFNKZcBRUUmi0%2FuhfuKGDKuATLMoMm2zDZS3wFku%2FmA%2FqiNaQnnvqRpVAdeZkzwoGPSuqeGkbRMGlF0pGxD7lD1MwnfykqQY6sQHvI3y4bBsZC0dj%2BxkOGD8b2qzxmjfpKP4%2Babs%2BnUyGVmTqQlYWYYvcAODWvwxLgN5%2BUbyhHR%2FRlETXABxyaf8DbiMcebo%2Bm8JaifeZfz653YOka99NOsO08iJafqPVeXAJ59OTXSeCugntpBJ4t53Ru5KlVpQXWMuPCYHz5zq9QWdz0Hie5nokW0r%2BYnAfSTIygQfYxIHUdIDydVCHSxehjYq6LEzJQQuTy3i%2FqudT1HE%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20231013T140326Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYW3DUAVV4%2F20231013%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=85f8397c1fe85b4cee4ae0fd07252d93dbf5b961fd0a805770600a0b3a05bdcc&hash=23d1ac222002fbdfd8694e459ed6ee7a97b1178d254705f0faa7acac2dc61f64&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=0043135472901820&tid=spdf-59375276-b0c2-44ef-8f81-57a2720b9927&sid=406f91a355a537440b0ac2509a9639992030gxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0f155a55505a505f515a07&rr=81582044bc21061e&cc=us
    ## Rxn kinematic_viscosity_water
    ##https://pubs.acs.org/doi/10.1021/acssuschemeng.2c03927
    ##https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.2009.0349


    ##Rxn Kinetics 

    Po = 350E-6 #atm
    K0(T,S) = exp(-60.2409 + 9345.17/T + 23.3585*log(0.01*T) + S*(0.023517-0.023656*0.01*T + 0.0047036*(0.01*T)^2)) #M/atm
    @show K0(292.5,25)
    K1(T,S) = 10^(-(3670.7/T-62.008+9.7944*log(T)-0.0118*S+0.000116*S^2)) #M
    K2(T,S) = 10^(-(1394.7/T + 4.777 - 0.0184*S + 0.000118*S^2)) #M
    KW = 1E-14 #M^2
    Kcal(T,S) = 10^(-(171.9065 + 0.077993*T - 2839.319/T - 71.595*log10(T) + (0.77712 - 0.0028426*T - 178.34/T)*S^(1/2) + 0.07711*S - 0.0041249*S^(3/2)))
    K_pos1(T) = exp(1246.98-(6.19*10^4)/T - 183.0*log(T))*3600 #hr-1
    K_neg1(T,S) = (K_pos1(T)/K1(T,S))/1E6 #hr-1 uM-1
    K_pos4(T,S) = (((499002.24*exp(4.2986E-04*S^2+5.75499E-05*S))*exp(-90166.83/(8.31451*T))/KW)*3600)/1E6 #uM-1 hr-1
    K_neg4(T,S) = K_pos4(T,S)*(KW/K1(T,S)) #hr -1
    r1(T,S,CO2) = -2*Kcal(T,S)/(CO2*(1/(molecular_weight_co2*1000))*K1(T,S)*K2(T,S))
    r2 = -1
    r3(T,S,CO2) = KW + CO2*K1(T,S)*(1/(molecular_weight_co2*1000)) #M^2
    r4(T,S,CO2) = 2*K1(T,S)*K2(T,S)*CO2*(1/(molecular_weight_co2*1000)) #M^3
   
    
    
    ##https://epic.awi.de/id/eprint/13960/1/Sch2006g.pdf
    ##file:///C:/Users/slant/Downloads/Dissociation_constants_of_carbonic_acid_in_seawate.pdf

    ## initial conditions for carbon species
    
    initial_pH = 1.00E-06 #M, initial pH = 8
    molecular_weight_h2co3 = 62.02/1000 #kg/mol
    molecular_weight_hco3 = 61.02/1000 #kg/mol
    molecular_weight_co3 = 60.01/1000 #kg/mol
    molecular_weight_ca = 40.08/1000 #kg/mol
    molecular_weight_caco3 = 100.9/1000 #kg/mol
    initial_h2co3 = 1.64E-06*1000*molecular_weight_h2co3 #kg/m3
    initial_hco3 = 3.28E-04*1000*molecular_weight_hco3 #kg/m3
    initial_co3 = 1.97E-08*1000*molecular_weight_co3 #kg/m3
    initial_ca = 400*(1/(1000*1000*molecular_weight_ca)) #mol/L
    initial_caco3 = 9.84E-09*1000*molecular_weight_caco3 #kg/m3


    co2_v = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8] #g/m3
    initial_co2_g(T,S) = co2_v[t]
    #Po*K0(T,S)*1000*1000*molecular_weight_co2 #g/m3
    typical_ow = 15.5E-6*1000*molecular_weight_co2 #kg/m3
    @show typical_ow
    ##https://www.iaea.org/sites/default/files/18/07/oa-chemistry-dickson-050916.pdf



    ## Mass Transfer Properties
    biomass_diffusion_coefficient_y = 1.0e-9 * 3600.0           #m^2/hour
    biomass_diffusion_coefficient_z = 1.0e-9 * 3600.0 * 100.0   #m^2/hour [100-fold higher in z-direction, due to added convection]

    ## Evaporation and Height Properties
    
    evaporation_constant(W) = 25 +19*W #kg dry air/m2-hr
    ambient_humidity_ratio(R_H, P) = 0.62198*(P*(R_H/100))/(reference_pressure-(P*(R_H/100))) #kg H2O/ kg dry air
    surface_humidity_ratio(T) = 0.62198*(saturated_vapor_pressure_water_air(T)/(reference_pressure-saturated_vapor_pressure_water_air(T))) #kg H2O/kg dry air
    evaporation_mass_flux(T,W,R_H,P) = -evaporation_constant(W)*(surface_humidity_ratio(T)-ambient_humidity_ratio(R_H, P)) #kg H2O/m2-hr
    evaporation_heat_flux(T,W,R_H,P) = evaporation_mass_flux(T,W,R_H,P)*heat_vaporization_water(T) #J/hr-m2
    ## above equations obtained from https://www.engineeringtoolbox.com/evaporation-water-surface-d_690.html 

    density_solution(T, S) = density_water(T) + S #kg/m3
    dVavgdx(T,S,R_H,W,P) = -((density_water(T)-density_water_vapor(T))*density_water(T)*specific_heat_capacity_water(T)*-evaporation_mass_flux(T,W,R_H,P))/(density_solution(T,S)*density_water_vapor(T)*density_solution(T,S)*specific_heat_capacity_solution(T, S)) #hr-1
    #change in Vavg with dx, derived from Wrobel, 2006
 
    height(x,T,S,W,R_H,P,H) = max(1/(((dynamic_viscosity_water(T)/density_solution(T,S))*(-evaporation_mass_flux(T,W,R_H,P))*(reactor_length/num_odes_y)*x)/((3*dynamic_viscosity_water(T)*volumetric_flow_rate_o)/(reactor_width*reactor_initial_liquid_level^3)) + (1/reactor_initial_liquid_level)), 0) #m
    #tot_energy_mod(H)*(reactor_length/proj_react(H))
    @show height(num_odes_y,298,24,3,44,saturated_vapor_pressure_water_air(298),0.25)
    
  

    ## Salinity Properties
    salinity_in = (1023.6 - density_water(298))
    salinity_o = (1023.6 - density_water(298))/density_water(298) #kg salt/kg water, obtained from "density of seawater @ 25 oC
    salinity(T,W,R_H,P,x,V) = salinity_o*density_water(298) + (salinity_o*density_water(298) - ((density_water(T)*(volumetric_flow_rate_o)*salinity_o)/(density_water(T)*(volumetric_flow_rate_o)-evaporation_mass_flux(T,W,R_H,P)*reactor_width*x*(reactor_length/num_odes_y)))*density_water(T)) #kg/m3



    params = (; num_odes_y,
                num_odes_z,
                time_end,
                time_interval,
                reference_pressure,
                gas_constant,
                stefan_boltzmann_constant,
                acceleration_gravity,
                molecular_weight_air,
                density_air,
                dynamic_viscosity_air,
                thermal_conductivity_air,
                specific_heat_capacity_air,
                prandtl_air,
                emissivity_air,
                molecular_weight_water,
                density_water,
                dynamic_viscosity_water,
                specific_heat_capacity_solution,                          
                specific_heat_capacity_water,
                thermal_conductivity_water,
                solar_reflectance_water,
                diffusion_coeff_water_air,
                saturated_vapor_pressure_water_air,
                heat_vaporization_water,
                emissivity_water,
                density_water_vapor,
                diffusion_coeff_co2_water,
                solubility_co2_water,
                molecular_weight_co2,
                thermal_diffusivity_concrete,
                thermal_conductivity_concrete,
                global_horizontal_irradiance_data,
                ambient_temperature_data,
                relative_humidity_data,
                wind_speed_data,
                data_begin,
                reactor_length,
                reactor_width,
                reactor_initial_liquid_level,
                reactor_depth,
                reactor_incline,
                rough_reynolds_number,
                boundary_layer_height_1,
                boundary_layer_height_2,
                boundary_layer_height_3,
                boundary_layer_height_4,
                velocity_profile_nw,
                velocity_profile_w,
                average_flow_velocity,
                velocity_profile_lam,
                average_flow_velocity_lam,
                volumetric_flow_rate,
                reynolds_number,
                input_hydraulic_diameter,
                input_temperature,
                ground_temperature,
                input_biomass_concentration,
                max_biomass_concentration,
                max_biomass_light_saturation,
                photosynthetic_efficiency,
                max_biomass_specific_growth_rate,
                salinity_factor,
                temperature_factor,
                pH_factor,
                Iave,
                phiL,
                co2_availability_factor,
                biomass_specific_growth_rate,
                co2_per_biomass,
                Po,
                K0,
                K1,
                K2,
                KW,
                Kcal,
                K_pos1,
                K_neg1,
                K_pos4,
                K_neg4,
                r1,
                r2,
                r3,
                r4,
                initial_pH,
                molecular_weight_h2co3,
                molecular_weight_hco3,
                molecular_weight_co3,
                molecular_weight_ca,
                molecular_weight_caco3,
                initial_h2co3,
                initial_hco3,
                initial_co3,
                initial_ca,
                initial_caco3,
                initial_co2_g,
                typical_ow,
                biomass_diffusion_coefficient_y,
                biomass_diffusion_coefficient_z,
                evaporation_constant,
                ambient_humidity_ratio,
                surface_humidity_ratio,
                evaporation_mass_flux,
                evaporation_heat_flux,
                density_solution,
                dVavgdx,
                height,
                salinity_in,
                salinity,
                )

    # Write parameters to file
    output = "name, value\n"
    for (name, val) in pairs(params)
        output = string(output, "$name, $val\n")
    end
    write("parameters_$filesuffix.csv", output)

    return params
end







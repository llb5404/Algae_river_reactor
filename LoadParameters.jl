
## LOAD PARAMETERS
## Options:
#   Load default parameters
#   Load parameters from an Excel file

using Debugger
using BasicInterpolators
break_on(:error)

function LoadDefaultParameters(filesuffix, l, q)
    ## PDE Discretization
    num_odes_y = 100
    num_odes_z = 11
    time_end = 732       #hours
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
    Ar(S) = 5.328 - 9.76 * 10 ^ (-2) * S + 4.04*10^(-4)*(S)^ 2 #unitless
    B(S) = -6.913 * 10 ^ (-3) + 7.351 * 10 ^ (-4) * (S) - 3.15*10^(-6)*(S)^2 #unitless
    C(S) = 9.6 * 10 ^ (-6) - 1.927 * 10 ^ (-6) * (S) + 8.23 * 10^(-9) *(S)^2 #unitless
    D(S) = 2.5 * 10 ^ (-9) + 1.666 * 10 ^ (-9) * (S) - 7.125 * 10^(-12)*(S)^2 #unitless
    specific_heat_capacity_solution(T,S) = 1000*(Ar(S) + B(S)*T + C(S)*(T^2) + D(S)*(T^3)) #J/kg/K
    #above information obtained from: 
    
    #K.G. Nayar, M.H. Sharqawy, L.D. Banchik, and J.H. Lienhard V, "Thermophysical properties of seawater: A review and new correlations that include pressure dependence," Desalination, Vol. 390, pp.1-24, 2016. doi:10.1016/j.desal.2016.02.024 (preprint)
    #Mostafa H. Sharqawy, John H. Lienhard V, and Syed M. Zubair, "Thermophysical properties of seawater: A review of existing correlations and data," Desalination and Water Treatment, Vol. 16, pp.354-380, April 2010. (PDF file which includes corrections through June 2017.)

    specific_heat_capacity_water(T) = specific_heat_capacity_solution(T,0) #J/kg/K             
    thermal_conductivity_water(T) = (0.565464 .+ 1.756E-03 .* T .- 6.46E-06 .* T.^2)*3600   # J/m/K/hr
    solar_reflectance_water = 0.1357 # fraction of light reflected
    diffusion_coeff_water_air(T) = 22.5E-06 .* ( (T .+ 273.15) ./ 273.15).^1.8*3600        # m^2 / sec
    saturated_vapor_pressure_water_air(T) = (exp(77.3450+0.0057*T-7235/T))/(T^8.2)  #Pascals
    heat_vaporization_water(T) = 2430159                                        # J / kg, T = 30 degC            [40660 * ( (647.3 - (T + 273.15)) / (647.3 - 373.2))^0.38 #J/mole, but isn't very accurate]
    emissivity_water = 0.97                                                     # dimensionless
    molecular_weight_water = 0.01801528                                         # kg/mole
    density_water_vapor(T) = 0.804 #kg/m3
    @show thermal_conductivity_water(298)
    # CO2 Properties
    diffusion_coeff_co2_water(T) = (13.942E-09*((293/227.0) - 1)^1.7094)*3600              # m^2 / hr [equation fitted based on data in https://pubs.acs.org/doi/10.1021/je401008s]
    solubility_co2_water(T) = (0.00025989 .* (T .- 273.15).^2 .- 0.03372247 .* (T .- 273.15) .+ 1.31249383) ./ 1000.0  #mole fraction of dissolved CO2 in water at equilibrium [equation fitted based on data from https://srd.nist.gov/JPCRD/jpcrd427.pdf]
    molecular_weight_co2 = 0.04401 #kg/mole                                             # kg / mole

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
    data_begin = 4345

    ## River Reactor Geometric Properties

    lengths = [25, 50, 75, 100, 125, 150, 175, 200]                           # meters
    reactor_length = 200
    real_length = lengths[q]

    reactor_width =  0.5                               # meters
    reactor_initial_liquid_level = 0.15             # meters, this is the sluice gate inlet height
    reactor_depth = 0   # meters
    height_diff = 0  # meters, 1 mm between start and end, practically 0

    evaporation_constant(W) = 25 +19*W #kg dry air/m2-hr
    ambient_humidity_ratio(R_H, P) = 0.62198*(P*(R_H/100))/(reference_pressure-(P*(R_H/100))) #kg H2O/ kg dry air
    surface_humidity_ratio(T) = 0.62198*(saturated_vapor_pressure_water_air(T)/(reference_pressure-saturated_vapor_pressure_water_air(T))) #kg H2O/kg dry air
    evaporation_mass_flux(T,W,R_H,P) = evaporation_constant(W)*(surface_humidity_ratio(T)-ambient_humidity_ratio(R_H, P)) #kg H2O/m2-hr
    #show for L = 50
    #flowrates from 25-200
    #L from 50 to 450 
    flow_rates = [5,10,15,20,25,30,35,40] #m3/hr
    length_flow = length(flow_rates)
    volumetric_flow_rate_o = flow_rates[l] #m^3/hr

    strip_flr = [8.5,7.5,6.5,5.5,4.5,3.5,2.5,1.5] #gpm (within range of low flr pumps on McMaster Carr)
    gpm = 0
    #strip_flr[t]
    volumetric_flow_rate_strip = gpm*0.00379*60 #m3/hr
    @show volumetric_flow_rate_strip
    strip_position = [trunc(Int,num_odes_y/2)] #halfway through reactor
    @show strip_position
    strip_concentration = 1500 #g/m3, saturated
    
    avg_velocity_o = volumetric_flow_rate_o/(reactor_initial_liquid_level*reactor_width)
    mass_o(T) = density_water(T)*reactor_initial_liquid_level*(reactor_length/num_odes_y)*reactor_width
    
    volumetric_flow_rate(T,W,R_H,P,x) = max((volumetric_flow_rate_o -evaporation_mass_flux(T,W,R_H,P)*(reactor_width/density_water(T))*x*(reactor_length/num_odes_y)),0)
    
    #+ volumetric_flow_rate_strip*(max(x-(strip_position[1]-1),0)/max(abs(x-(strip_position[1]-1)),1)) ) 

    term(V,H) = V/((2/3)*H^2)
    velocity_profile_lam(V,y,H) = term(V,H)*((H-y*(H/num_odes_z))*H + 0.5*(H-y*(H/num_odes_z))^2)
    hydraulic_diameter(H) = 4 * (reactor_width * H) / (reactor_width + 2 * H) #m
    reynolds_number(H,T,V) = V .* density_water(T) .* hydraulic_diameter(H) ./ (dynamic_viscosity_water(T))
    
    
    input_hydraulic_diameter = 4 * (reactor_width * reactor_initial_liquid_level) / (reactor_width + 2 * reactor_initial_liquid_level) # meters
    input_temperature = 293.15                          # Kelvin
    ground_temperature = 290.0                          # Kelvin

    ## Biomass Properties
    input_biomass_concentration = 3000 / 1000.0         # kg/m^3
    max_biomass_concentration = 10000.0 / 1000.0        # kg/m^3 where light is 100% absorbed
    max_biomass_light_saturation = 900.0 / 4.57         # 900 μmol m−2 s−1 converted to W/m^2
    max_biomass_specific_growth_rate = 7.9/24   # 1 / hour (obtained from Krishnan at T_opt, 35 salinity)
    
    #biomass growth adjustment factors
    salinity_factor(S) = -0.0002*S^2 + 0.0107*S + 0.8661 #unitless, S in kg/kg sw
    ## obtained from 2-degree polynomial fit of data obtained from Krishnan et al.
    T_opt = 308 #in C, optimum temperature
    T_min = 277.2 #in C, minimum temperature
    T_max = 320 #in C, maximum temperature
    temperature_factor_g(T) = (T_opt-T_min)*(T-T_opt) #unitless
    temperature_factor_f(T) = (T_opt-T_max)*(T_opt+T_min-2*(T)) #unitless
    temperature_factor(T) = ((T-T_max)*(T-T_min)^2)/((T_opt-T_min)*(temperature_factor_g(T)-temperature_factor_f(T))) #unitless
    


    ## opt, min, max temperatures given by Krishnan et al 
    ##temperature factor obtained from model given by Greene et al https://www.osti.gov/servlets/purl/1806213

    photosynthetic_efficiency = 0.025                   # fraction of sunlight converted to chemical energy during photosynthesis
    Vol(dz) = (reactor_length/num_odes_y)*dz*reactor_width #m3
    
    co2_availability_factor(C_co2) = 1
    #((130.43*C_co2)/(4286.4+C_co2+(C_co2^2/0.012898)))*(1/0.114)


    
    #(C_co2) ./ ((C_co2) + 4.26 + (C_co2)^2/(250)) #unitless, input is CO2 in g/m3

    pH_opt = 7.0

    pH_CO2(pH) = exp((7.7291-pH + (pH_opt-6.328))/0.402)
    pH_factor(pH) = pH_CO2(pH)/(pH_CO2(pH) + 4.26 + (pH_CO2(pH)^2)/250)

    biomass_specific_growth_rate(T, S, C_co2) = max_biomass_specific_growth_rate .* co2_availability_factor(C_co2).*temperature_factor(T).*salinity_factor(S) #1/hr
    
    
    co2_per_biomass = 1.88   #kg CO2/kg algae

    co2_v = [0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10]

    sparge_fpm = 12.5 #fpm, based on standard sparger design range https://mottcorp.com/wp-content/uploads/2020/05/Sparger-Design-Guide.pdf
    
    P_atm = 1 #atm
    PO3 = (0.1E-03)/1000 #mol/kg soln, https://resourcewatch.org/data/explore/f1aa9ec7-c3b6-441c-b395-96fc796b7612?section=Discover&selectedCollection=&zoom=2.422253880286214&lat=51.07099144291875&lng=-85.84319789585153&pitch=0&bearing=0&basemap=dark&labels=light&layers=%255B%257B%2522dataset%2522%253A%2522f1aa9ec7-c3b6-441c-b395-96fc796b7612%2522%252C%2522opacity%2522%253A1%252C%2522layer%2522%253A%25221122cdbf-cb73-467a-bb25-ad86ac491136%2522%257D%255D&aoi=&page=1&sort=most-viewed&sortDirection=-1
    Si = 0.78E-06 #mol/kg soln, https://plymsea.ac.uk/id/eprint/1451/1/The_determination_of_silicate_in_sea_water.pdf
    NH4 = 30E-09

    co2_to_M = (1/(molecular_weight_co2*1000*1000))

    pH_inter_data = zeros(1001)

    csv_reader1 = CSV.File("pH_inter.csv", header=["col1"])
    for (i,row) in enumerate(csv_reader1)
        pH_inter_data[i] = convert(Float64, row.col1)   
     
    end
 
    vector1 = range(1E-09,70001,1001) #uses DIC

    pH_interp = LinearInterpolator(vector1, pH_inter_data) #function that uses DIC as input, pH as output

    CO2_inter_data = zeros(1001)
    csv_reader2 = CSV.File("CO2_inter.csv", header=["col1"])

    for (i,row) in enumerate(csv_reader2)
        CO2_inter_data[i] = convert(Float64, row.col1)    #g/m3
    end

    CO2_interp = LinearInterpolator(vector1, CO2_inter_data) #function that uses DIC as input, CO2 in mol/L as output
    
    ## Mass Transfer Properties
    biomass_diffusion_coefficient_y = 1.0e-9 * 3600.0           #m^2/hour
    biomass_diffusion_coefficient_z = 1.0e-9 * 3600.0           #m^2/hour
    
    ## Evaporation and Height Properties
    
    evaporation_heat_flux(T,W,R_H,P) = evaporation_constant(W)*(surface_humidity_ratio(T)-ambient_humidity_ratio(R_H, P))*heat_vaporization_water(T) #J/hr-m2
    ## above equations obtained from https://www.engineeringtoolbox.com/evaporation-water-surface-d_690.html 

    density_solution(T, S) = density_water(T) + 0.765*S #kg/m3
    dVavgdx(T,S,R_H,W,P) = ((density_water(T)-density_water_vapor(T))*density_water(T)*specific_heat_capacity_water(T)*evaporation_mass_flux(T,W,R_H,P))/(density_solution(T,S)*density_water_vapor(T)*density_solution(T,S)*specific_heat_capacity_solution(T, S)) #hr-1
    avg_velocity(T,S,R_H,W,P,x) = max(avg_velocity_o + dVavgdx(T,S,R_H,W,P)*x*(reactor_length/num_odes_y), 1E-09)

    height(T,S,R_H,W,P,x) = volumetric_flow_rate(T,W,R_H,P,x)/(S*reactor_width)

    ## Salinity Properties
    salinity_in = 35 #g/kg sw

    ##co2 properties
    mol_frac_co2 = 0.0 #mol frac of CO2 in gas phase
    mol_frac_air = 1-mol_frac_co2
    kg_mol_gas = mol_frac_co2*(1/molecular_weight_co2) + mol_frac_air*(1/molecular_weight_air) #kg gas/mol
    floor_oc_coeff = 0
    #co2_v[t] #fraction of floor area occupied by spargers
    tot_area = reactor_width*(reactor_length)*floor_oc_coeff
    mass_flr_gas = (sparge_fpm*(60/3.28)*tot_area)*density_air(292.5) #kg gas/hr
    G = (mass_flr_gas/kg_mol_gas) #molar flow of gas, mol/hr
    henry_const(T) = 0.035*exp(2400*((1/T) - (1/298.5)))*density_water(T)*(1/0.987)*(1/1000) #henry' constant (mol/L/atm)
    kla_CO2 = 0.0020*3600 #hr-1
    A = tot_area/num_odes_y #m2
    y_out(T,C,H,mf) = (1/henry_const(T))*(C*co2_to_M+(mf*henry_const(T) - C*co2_to_M)*exp(-kla_CO2*(A/G)*H*henry_const(T)*1000))
    dMt(T,C,H) = G*(mf - y_out(T,C,H,mf))
    #https://www.sciencedirect.com/science/article/pii/S0960852412012047?casa_token=Dg_MAh0F[%E2%80%A6]RwNik8I_cva5L1jX7aB20_ytLrzqUqHu6U7HcAf2xvkFgdibxBymq8QiTI

    co2_init = 22.5 #g/m3
    DIC_init = 2002+co2_init*(1/(44.01*1000))*1E06 #uM
    pH_init = pH_interp(DIC_init)
    @show pH_init
    ratio_init = (co2_init*(1/molecular_weight_co2))/DIC_init

    salinity_o = salinity_in #kg salt/kg water, obtained from "density of seawater @ 25 oC
    #salinity(T,W,R_H,P,x,V) = density_water(T)*salinity_o*(((density_water(T)*(volumetric_flow_rate_o))/(density_water(T)*(volumetric_flow_rate_o)-evaporation_mass_flux(T,W,R_H,P)*reactor_width*x*(reactor_length/num_odes_y)))) #kg/m3
    
    salinity(T,R_H,W,P,x) = (salinity_o*volumetric_flow_rate_o)/(volumetric_flow_rate(T,W,R_H,P,x)) #kg salt in a section/kg sw


    karmen_c = 0.41
    shear_velocity(T,S,R_H,W,P,x) = max(avg_velocity(T,S,R_H,W,P,x)/(2.5*log((12.14*height(T,S,R_H,W,P,x))/karmen_c)),0)
    velocity_profile_turb(T,S,R_H,W,P,x,y,H) = max(shear_velocity(T,S,R_H,W,P,x)*2.5*log((33*(H - y*(H/num_odes_z)))/karmen_c),0)
    
    
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
                lengths,
                reactor_length,
                reactor_width,
                reactor_initial_liquid_level,
                reactor_depth,
                height_diff,
                evaporation_constant,
                ambient_humidity_ratio,
                surface_humidity_ratio,
                evaporation_mass_flux,
                flow_rates,
                length_flow,
                avg_velocity_o,
                mass_o,
                avg_velocity,
                volumetric_flow_rate,
                height,
                velocity_profile_lam,
                hydraulic_diameter,
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
                co2_availability_factor,
                biomass_specific_growth_rate,
                co2_per_biomass,
                co2_v,
                DIC_init,
                ratio_init,
                P_atm,
                PO3,
                Si,
                NH4,
                pH_interp,
                CO2_interp,
                biomass_diffusion_coefficient_y,
                biomass_diffusion_coefficient_z,
                evaporation_heat_flux,
                density_solution,
                G,
                y_out,
                dMt,
                salinity_in,
                mol_frac_co2,
                co2_init,
                salinity,
                dVavgdx,
                velocity_profile_turb,
                pH_init,
                volumetric_flow_rate_o,
                strip_position,
                volumetric_flow_rate_strip,
                strip_concentration,
                strip_flr,
                real_length
                )

    # Write parameters to file
    output = "name, value\n"
    for (name, val) in pairs(params)
        output = string(output, "$name, $val\n")
    end
    write("parameters_$filesuffix.csv", output)

    return params
end







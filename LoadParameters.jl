cd("C:/Users/slant/OneDrive/Desktop/Julia_2")
## LOAD PARAMETERS
## Options:
#   Load default parameters
#   Load parameters from an Excel file  <to-do>

using Debugger

break_on(:error)

function LoadDefaultParameters(filesuffix)
    ## PDE Discretization
    num_odes_y = 50
    num_odes_z = 5
    time_end = 960.0        #hours
    time_interval = 1.0     #hours
    ## Physical Constants
    reference_temperature = 20.0                                        # deg Celsius
    reference_pressure = 101315.0                                       # Pascals
    gas_constant = 8.31446261815324                                     # J / K / mole
    stefan_boltzmann_constant = 5.67037442E-08                          # W / m^2 / K^4
    acceleration_gravity = 9.81                                         # m/s2
    # Air Properties
    density_air(T) = reference_pressure ./ (287.058 .* T)                 # kg/m3
    dynamic_viscosity_air(T) = 1.458E-06 .* (T.^1.5) ./ (T .+ 110.4)        # Pa sec [kg/m/sec]
    thermal_conductivity_air(T) = 0.02624 .* (T/300).^0.8646              # W/m*K
    specific_heat_capacity_air(T) = 1002.5 .+ (275E-06) .* (T .- 200).^2      # J/kg*K
    prandtl_air(T) = specific_heat_capacity_air(T) .* dynamic_viscosity_air(T) ./ thermal_conductivity_air(T)             # dimensionless
    emissivity_air = 0.80                                               # dimensionless
    molecular_weight_air = 0.02896                                      # kg/mole [average of air mixture]

    # Water Properties
    density_water_20degC = 998.2071                                             # kg / m^3
    expansion_coefficient_water = 0.0002                                        # 1/degC
    density_water(T) = density_water_20degC ./ (1 .+ expansion_coefficient_water .* (T .- reference_temperature))   #kg/m^3
    dynamic_viscosity_water(T) = 2.414E−05 .* 10^( 247.8 ./ (T .+ 133.15))         # Pa sec [kg/m/sec]
    specific_heat_capacity_water(T) = 4184.0                                    # J/kg/K
    thermal_conductivity_water(T) = 0.565464 .+ 1.756E-03 .* T .- 6.46E-06 .* T.^2   # J/m/K/sec
    solar_reflectance_water = 0.1357                                            # fraction of light reflected
    diffusion_coeff_water_air(T) = 22.5E-06 .* ( (T .+ 273.15) ./ 273.15).^1.8        # m^2 / sec
    saturated_vapor_pressure_water_air(T) = (exp(77.3450+0.0057*T-7235/T))/(T^8.2)  #Pascals
    heat_vaporization_water(T) = 2430159                                        # J / kg, T = 30 degC            [40660 * ( (647.3 - (T + 273.15)) / (647.3 - 373.2))^0.38 #J/mole, but isn't very accurate]
    emissivity_water = 0.97                                                     # dimensionless
    molecular_weight_water = 0.01801528                                         # kg/mole
    
    # CO2 Properties
    diffusion_coeff_co2_water(T) = (0.0816 .* T .- 22.601) .* 1.0E-09              # m^2 / sec [equation fitted based on data in https://pubs.acs.org/doi/10.1021/je401008s]
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
    reactor_length = 200.0                              # meters
    reactor_width =  20.0                               # meters
    reactor_initial_liquid_level = 2.0                 # meters
    reactor_depth = 2E-06 * reactor_initial_liquid_level    # meters
    reactor_incline = (asin(reactor_depth/reactor_length))  # radians
    @show reactor_incline
    ## River Reactor Operating Parameters
    input_volumetric_flow_rate = 5000     # m^3/hour
    input_max_flow_velocity = 2.0 * input_volumetric_flow_rate / reactor_width / reactor_initial_liquid_level   # m/hour
    input_average_flow_velocity = ((density_water(20)*acceleration_gravity*tan(reactor_incline)*reactor_initial_liquid_level^2)/(3*dynamic_viscosity_water(20)))*3600 #m/hr
    change_potential_energy = -2*sin(reactor_incline)*(reactor_length/num_odes_y)*acceleration_gravity*(3600^2) #in m2/hr2
    initial_average_flow_velocity(x) = sqrt((input_average_flow_velocity)^2+(change_potential_energy*x)) #average_flow_velocity constant despite height change
    @show initial_average_flow_velocity(0)
    initial_height_profile(x) = (input_average_flow_velocity*reactor_initial_liquid_level)/initial_average_flow_velocity(x)
    ## not initial
    velocity_profile(z, H, T) = (density_water(T)*acceleration_gravity*tan(reactor_incline)*z*(2*H-z))/(2*dynamic_viscosity_water(T))*3600 #z = dz(i)*j
    average_flow_velocity(H, T) = ((density_water(T)*acceleration_gravity*tan(reactor_incline)*H^2)/(3*dynamic_viscosity_water(T)))*3600
    
    
      
    input_hydraulic_diameter = 4 * (reactor_width * reactor_initial_liquid_level) / (reactor_width + 2 * reactor_initial_liquid_level) # meters
    input_temperature = 293.15                          # Kelvin
    ground_temperature = 290.0                          # Kelvin

    #Reynolds number for a duct (open surface), Re = Vavg * density * hydraulic diameter / viscosity
    input_reynolds_number(T) = input_average_flow_velocity .* density_water(T) .* input_hydraulic_diameter ./ viscosity_water(T)   # dimensionless

    ## Biomass Properties
    input_biomass_concentration = 10.0 / 1000.0         # kg/m^3
    max_biomass_concentration = 10000.0 / 1000.0        # kg/m^3 where light is 100% absorbed
    max_biomass_light_saturation = 900.0 / 4.57         # 900 μmol m−2 s−1 converted to W/m^2
    max_biomass_specific_growth_rate = log(2.0) / 3.0   # 1 / hour
    photosynthetic_efficiency = 0.025                   # fraction of sunlight converted to chemical energy during photosynthesis
    threshold_dissolved_co2_growth = 0.5                # kg CO2 / m^3 water, minimum dissolved CO2 concentration before growth begins to slow
    
    Vol(dz) = (reactor_length/num_odes_y)*dz*reactor_width
    Iave(M_biomass_z, GHI, z, dz) = GHI .* 0.45 .* (1 - min(sum(M_biomass_z[1:z]/(Vol(dz))).*dz./max_biomass_concentration,1.0) ) #
    phiL(M_biomass_z, GHI, z, dz) = Iave(M_biomass_z, GHI, z, dz) .* exp(1 - Iave(M_biomass_z, GHI, z, dz)./max_biomass_light_saturation) ./ max_biomass_light_saturation
    co2_availability_factor(C_co2, dz) = min.(C_co2/Vol(dz), threshold_dissolved_co2_growth) ./ threshold_dissolved_co2_growth
    biomass_specific_growth_rate(M_biomass_z, GHI, z, dz, C_co2) = max_biomass_specific_growth_rate .* phiL(M_biomass_z, GHI, z, dz) .* co2_availability_factor(C_co2,dz)
    co2_per_biomass = 0.70   #based on 2.661 kg CO2 emitted when burning 1 gallon of algae (about 3.79 kg)

    ## Mass Transfer Properties
    biomass_diffusion_coefficient_y = 1.0e-9 * 3600.0           #m^2/hour
    biomass_diffusion_coefficient_z = 1.0e-9 * 3600.0 * 100.0   #m^2/hour [100-fold higher in z-direction, due to added convection]

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
                specific_heat_capacity_water,
                thermal_conductivity_water,
                solar_reflectance_water,
                diffusion_coeff_water_air,
                saturated_vapor_pressure_water_air,
                heat_vaporization_water,
                emissivity_water,
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
                reactor_incline,
                reactor_depth,
                input_volumetric_flow_rate,
                input_max_flow_velocity,
                input_average_flow_velocity,
                change_potential_energy,
                initial_average_flow_velocity,
                initial_height_profile,
                velocity_profile,
                average_flow_velocity,
                input_hydraulic_diameter,
                input_reynolds_number,
                input_temperature,
                ground_temperature,
                input_biomass_concentration,
                max_biomass_concentration,
                max_biomass_light_saturation,
                photosynthetic_efficiency,
                max_biomass_specific_growth_rate,
                Vol,
                Iave,
                phiL,
                co2_availability_factor,
                biomass_specific_growth_rate,
                co2_per_biomass,
                biomass_diffusion_coefficient_y,
                biomass_diffusion_coefficient_z,
                )

    # Write parameters to file
    output = "name, value\n"
    for (name, val) in pairs(params)
        output = string(output, "$name, $val\n")
    end
    write("parameters_$filesuffix.csv", output)

    return params
end







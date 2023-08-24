cd("C:/Users/slant/OneDrive/Desktop/Julia_2")
function HeightChange!(dX, H, Temperature, params, t)
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data
    data_begin = params.data_begin

    t_hour1 = floor(Int64, t)
    t_hour2 = floor(Int64, t)+1

    GHI = GHI_Data[data_begin + t_hour1] * (t-t_hour1) + GHI_Data[data_begin + t_hour2] * (t_hour2 - t)
    WNDSPD = max.((WNDSPD_Data[data_begin + t_hour1] * (t-t_hour1) + WNDSPD_Data[data_begin + t_hour2] * (t_hour2 - t)),0)
    RH = RH_Data[data_begin + t_hour1] * (t-t_hour1) + RH_Data[data_begin + t_hour2] * (t_hour2 - t)
    Tamb = max.(((Tamb_Data[data_begin + t_hour1] * (t-t_hour1) + Tamb_Data[data_begin + t_hour2] * (t_hour2 - t))),0)

    rho_water = params.density_water(Tamb)
    rho_air = params.density_air(Tamb)
    rho_sol = params.density_solution
    dyn_visc_water = params.dynamic_viscosity_water(Tamb)
    dyn_visc_air = params.dynamic_viscosity_air(Tamb)
    cp_water = params.specific_heat_capacity_water(Tamb)
    cp_air = params.specific_heat_capacity_air(Tamb)
    pr_air = params.prandtl_air(Tamb)
    P_w = params.saturated_vapor_pressure_water_air
    P_a = params.saturated_vapor_pressure_water_air(Tamb)
    k_water = params.thermal_conductivity_water(Tamb)
    k_air = params.thermal_conductivity_air(Tamb)
    D_w_a = params.diffusion_coeff_water_air(Tamb)
    hfg_water = params.heat_vaporization_water(Tamb)
    g = params.acceleration_gravity
    theta = params.reactor_incline

    alpha_w = params.solar_reflectance_water
    M_water = params.molecular_weight_water
    R = params.gas_constant
    sigma = params.stefan_boltzmann_constant
    epsilon_water = params.emissivity_water
    epsilon_air = params.emissivity_air
    diff_concrete = params.thermal_diffusivity_concrete
    ground_temperature = params.ground_temperature
    k_concrete = params.thermal_conductivity_concrete
    f_a = params.photosynthetic_efficiency

    Q = params.volumetric_flow_rate       # m^3/hour
    L = params.reactor_length                   # m
    W = params.reactor_width                    # m
    Vavg(H) = params.average_flow_velocity(H,Tamb) #m/hour
   
    P_atm = params.reference_pressure
    Ny = params.num_odes_y
    Nz = params.num_odes_z

    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    

    Nelements = (Ny+1) * (Nz+1)
    H = max.(H, 0.0) #don't allow negative values
    dH = zeros( Nelements, 1)


    K(T,S) = params.dVavgdx(T,S,RH,WNDSPD,P_a)
    M_Evap(T) = params.evaporation_mass_flux(T,WNDSPD,RH,P_a)

    
  

    for i = 1:Ny
       
 
        # BC 2
        dH[pos2idx(i, 0)] = (((dyn_visc_water/rho_water)*(M_Evap(Temperature[pos2idx(i,0)])+K(Temperature[pos2idx(i,0)],0)*rho_sol(Temperature[pos2idx(i,0)],0)))/(rho_water*g*tan(theta)*H[pos2idx(i,0)]^2))*Vavg(H[pos2idx(i,0)])
        # BC 3
        dH[pos2idx(i, Nz)] = ((dyn_visc_water/rho_water)*(M_Evap(Temperature[pos2idx(i,0)])+K(Temperature[pos2idx(i,Nz)],0)*rho_sol(Temperature[pos2idx(i,Nz)],0))/(rho_water*g*tan(theta)*H[pos2idx(i,Nz)]^2))*Vavg(H[pos2idx(i,Nz)])
        

    end

    for j = 0:Nz

        dH[pos2idx(0, j)] = 0.0 #((k_water * (T[pos2idx(0, max(0,j-1))] - 2*T[pos2idx(0, j)] + T[pos2idx(Ny, min(Nz,j+1))]) / dz^2))/(rho_water*cp_water)
    end

    for i=1:Ny
        for j=1:Nz-1
            dH[pos2idx(i,j)] = ((dyn_visc_water/rho_water)*(M_Evap(Temperature[pos2idx(i,0)])+K(Temperature[pos2idx(i,j)],0)*rho_sol(Temperature[pos2idx(i,j)],0))/(rho_water*g*tan(theta)*H[pos2idx(i,j)]^2))* Vavg(H[pos2idx(i,j)])
        
        end
    end

    @views dX[3*Nelements+1:4*Nelements] .= dH        # order matters!  The @views operator takes a slice out of an array without making a copy.
    nothing
end

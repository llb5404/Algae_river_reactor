
function HeatTransfer!(dX, T, Height, params, t)
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
    rho_solution(S) = params.density_solution(Tamb,S)
    rho_air = params.density_air(Tamb)
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
    V_profile(z,H) = params.velocity_profile(z,H,Tamb)
    P_atm = params.reference_pressure
    Ny = params.num_odes_y
    Nz = params.num_odes_z

    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    dy = L / Ny

    Nelements = (Ny+1) * (Nz+1)
    T = max.(T, 0.0) #don't allow negative values
    dT = zeros( Nelements, 1)
    H_o = params.reactor_initial_liquid_level
    K(T,S) = params.dVavgdx(T,S,RH,WNDSPD,P_a)
    M_Evap(T) = params.evaporation_mass_flux(T,WNDSPD,RH,P_a)

    Hght(x,T,S) = params.height(x,T,S,WNDSPD,RH,P_a,H_o)
    dz(x,T,S) = Hght(x,T,S)/Nz
   
    

    #constants in all heat flux functions
    #This function defines the heat flux from direct solar radiation

    Q_Solar = (1-0.03) * alpha_w * GHI * 3600;                                  #goes in, positive. removed * (1-f_a) as the fraction of biomass converting light into chemical energy is very small
    #Q_Solar = (1-0.03) * alpha_w * (1-f_a) * GHI * 3600;                       #goes in, positive

    # Source: Yadala and Cremaschi, 2016
    #This function defines the heat flux from reradiation from the pond (Q_out)

    Q_Rerad(Temp) = -1.0 * sigma * epsilon_water .* Temp.^4 * 3600              #goes out, negative

    #This function defines the heat flux from longwave atmospheric radiation
    Q_Longwave_Atmo = epsilon_water * epsilon_air * sigma .* Tamb.^4 * 3600     #goes in, positive

    #This function calculates the evaporation rate of the pond as well as the
    #cooling effect due to that evaporation. Straight from Yadala and Cremaschi, 2016.

    #Calculate the Reynold's Number
    kin_visc_a = dyn_visc_air / rho_air      # m^2/s
    L_c = L
    Re_L = L_c * WNDSPD / kin_visc_a

    Q_Evap(Temp) = params.evaporation_heat_flux(Temp,WNDSPD,RH,P_a)

    #Q_Recover(Temp) = M_Evap(Temp) * cp_water * (293.15-Temp) * 3600       #not sure what this is

    #This function calculates the evaporation rate of the pond as well as the
    #cooling effect due to that evaporation. Straight from Yadala and
    #Cremaschi, 2016.

    #This function calculates the convection heat flux within the pond

    Nu_L = 0.035 * (Re_L^0.8) * pr_air^(1/3)

    #Caclculate the convection coefficient given the Nusselt number
    h_conv = Nu_L * k_air / L_c

    #Calculate the convective heat transfer
    Q_Conv(Temp) = h_conv * (Tamb - Temp) * 3600                                              # J/hr/m2     goes in if Tamb>Ty - positive

    # Pretty much straight from Yadala and Cremaschi, 2016. Except for the added
    # correlation for laminar flow and averaging the two if in the transisition period.

    # Heat Transfer through Ground
    l_ref = 4400 * diff_concrete^0.5
    Q_Ground(Temp) = -k_concrete * (Temp - ground_temperature) / l_ref * 3600                 #J/hr/m2     goes out if Temp > ground temperature

    Q_sum1(Temp) = Q_Longwave_Atmo + Q_Rerad(Temp) + Q_Conv(Temp) + Q_Evap(Temp) + Q_Solar
    Q_sum2(Temp) = Q_Ground(Temp)
   
    #BC1: T(y,z) at z = 0 is constant [Tin]
    #BC2: dT(y,z) at z = 0 includes the d2T/dy2, Vy*dT/dy, and Q terms (Qsol, Qrerad, Qlw, Qcov, Qevap)
    #BC3: dT(y,z) at z = H includes the d2T/dy2, Vy*dT/dy, and Q terms (Qground)
    # rho_water*cp_water*dT/dt = k*[ d2T/dy2 + d2T/dz2 ] - Vy * dT/dy + Q
    #0 to Ny?

    ##T(-1,j) = T(0,j)
    ##T(Ny+1,j) = T(Ny,j)

    
  

    for i = 1:Ny-1
       
 
        # BC 2
        dT[pos2idx(i, 0)] = (  k_water*(T[pos2idx(i-1, 0)] - 2*T[pos2idx(i, 0)] + T[pos2idx(min(Ny,i+1), 0)]) / dy^2    #conduction in y-dir
                             + k_water*(T[pos2idx(i, 0)] + T[pos2idx(i,0+1)] - 2*T[pos2idx(i, 0)]) / dz(i,T[pos2idx(i, 0)],0)^2               #conduction in z-dir at boundary
                             + Q_sum1(T[pos2idx(i,0)]) * dy * W                                                         #heat generation at boundary
                             - V_profile(0, Hght(i,T[pos2idx(i, 0)],0)) * ( T[pos2idx(i,0)] - T[pos2idx(i-1,0)]) / dy                                       #heat convection at boundary
                             ) / (rho_water * cp_water)

        # BC 3
        dT[pos2idx(i, Nz)] = (  k_water*(T[pos2idx(i-1, Nz)] - 2*T[pos2idx(i, Nz)] + T[pos2idx(min(Ny,i+1), Nz)]) / dy^2 #conduction in y-dir
                              + k_water*(T[pos2idx(i, Nz-1)] + T[pos2idx(i,Nz)] - 2*T[pos2idx(i, Nz)]) / dz(i,T[pos2idx(i, Nz)],0)^2            #conduction in z-dir at boundary
                              + Q_sum2(T[pos2idx(i,Nz)]) * dy * W                                                        #heat generation at boundary
                              - V_profile(dz(i,T[pos2idx(i, Nz)],0)*Nz, Hght(i,T[pos2idx(i, Nz)],0)) * ( T[pos2idx(i,Nz)] - T[pos2idx(i-1,Nz)]) / dy                                        #heat convection at boundary
                             ) / (rho_water*cp_water)
    end

    for j = 0:Nz

        dT[pos2idx(0, j)] = 0 #((k_water * (T[pos2idx(0, max(0,j-1))] - 2*T[pos2idx(0, j)] + T[pos2idx(Ny, min(Nz,j+1))]) / dz^2))/(rho_water*cp_water)

        dT[pos2idx(Ny, j)] = (  k_water * (T[pos2idx(Ny, max(0,j-1))] - 2*T[pos2idx(Ny,j)] + T[pos2idx(Ny, min(Nz,j+1))]) / dz(Ny,T[pos2idx(Ny, j)],0)^2    #conduction in z-dir
                              + k_water * (T[pos2idx(Ny-1, j)] - 2*T[pos2idx(Ny,j)] + T[pos2idx(Ny,j)]) / dy^2                      #conduction in y-dir at boundary
                              - V_profile(dz(Ny,T[pos2idx(Ny, j)],0)*j, Hght(Ny,T[pos2idx(Ny, j)],0)) * (T[pos2idx(Ny,j)] - T[pos2idx(Ny-1,j)]) / dy                                      #heat convection at boundary
                             ) / (rho_water*cp_water)
    end

    for i=1:Ny-1
        for j=1:Nz-1
            dT[pos2idx(i,j)] = (   k_water * (T[pos2idx(i-1, j)] - 2*T[pos2idx(i, j)] + T[pos2idx(i+1, j)]) / dy^2      #conduction in y-dir
                                 + k_water * (T[pos2idx(i, j-1)] - 2*T[pos2idx(i, j)] + T[pos2idx(i, j+1)]) / dz(i,T[pos2idx(i, j)],0)^2      #conduction in z-dir
                                 - V_profile(dz(i,T[pos2idx(i, j)],0)*j, Hght(i,T[pos2idx(i, j)],0)) * (T[pos2idx(i,j)] - T[pos2idx(i-1,j)]) / dy                         #heat convection
                               ) / (rho_water*cp_water)
        end
    end


    #@show dT
    @views dX[Nelements+1:2*Nelements] .= dT        # order matters!  The @views operator takes a slice out of an array without making a copy.
    nothing
end
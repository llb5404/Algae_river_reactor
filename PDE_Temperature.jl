
function HeatTransfer!(dX, T, params, t)
    #Units for params can be found in LoadParameters.jl
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
    dyn_visc_air = params.dynamic_viscosity_air(Tamb)
    cp_water = params.specific_heat_capacity_water(Tamb)
    pr_air = params.prandtl_air(Tamb)
    P_a = params.saturated_vapor_pressure_water_air(Tamb)
    k_water = params.thermal_conductivity_water(Tamb)
    k_air = params.thermal_conductivity_air(Tamb)

    alpha_w = params.solar_reflectance_water
    sigma = params.stefan_boltzmann_constant
    epsilon_water = params.emissivity_water
    epsilon_air = params.emissivity_air
    diff_concrete = params.thermal_diffusivity_concrete
    ground_temperature = params.ground_temperature
    k_concrete = params.thermal_conductivity_concrete

    L = params.reactor_length                   
    W = params.reactor_width                    
    H_o = params.reactor_initial_liquid_level
    Ny = params.num_odes_y
    Nz = params.num_odes_z

    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    dy = L / Ny

    Nelements = (Ny+1) * (Nz+1)
    T = max.(T, 0.0) #don't allow negative values
    dT = zeros( Nelements, 1)
    K(T,S) = params.dVavgdx(T,S,RH,WNDSPD,P_a)
    M_Evap(T) = params.evaporation_mass_flux(T,WNDSPD,RH,P_a)

    Hght(x,T,S) = params.height(x,T,S,WNDSPD,RH,P_a,H_o)
    dz(x,T,S) = Hght(x,T,S)/Nz

   ## Velocity Profile

    Re(H,T,V) = params.reynolds_number(H,T,V)
    Re_star(H,T) = params.rough_reynolds_number(H,T)
    V_profile = zeros(Nelements, 1)
    Vavg_lam(H,T) = params.average_flow_velocity_lam(H,T)
    Vavg(H,B) = params.average_flow_velocity(H,B)
    S(T,x,V) = params.salinity(T,WNDSPD,RH,P_a,x,V)
    dz_v = zeros(Nelements,1)
    Sal = zeros(Nelements, 1)
    Ht = zeros(Nelements,1)

    for i in 0:Ny
        for j in 0:Nz
            #cut-off for laminar flow = 500
            if Re(Hght(i,T[pos2idx(i,j)],S(T[pos2idx(i,j)],i,V_profile[pos2idx(i,j)])),T[pos2idx(i,j)],Vavg_lam(Ht[pos2idx(i,j)],T[pos2idx(i,j)])) > 500
                #Turbulent Flow: choose boundary layer thickness based on "rough reynolds number"
                if Re_star(Hght(i,T[pos2idx(i,j)],S(T[pos2idx(i,j)],i,V_profile[pos2idx(i,j)])),T[pos2idx(i,j)]) <= 4
                    global Bd = params.boundary_layer_height_1(Ht[pos2idx(i,j)],T[pos2idx(i,j)]) #boundary layer, m
                elseif Re_star(Hght(i,T[pos2idx(i,j)],S(T[pos2idx(i,j)],i,V_profile[pos2idx(i,j)])),T[pos2idx(i,j)]) <= 11
                    global Bd = params.boundary_layer_height_2(Ht[pos2idx(i,j)],T[pos2idx(i,j)]) #boundary layer, m
                elseif Re_star(Hght(i,T[pos2idx(i,j)],S(T[pos2idx(i,j)],i,V_profile[pos2idx(i,j)])),T[pos2idx(i,j)]) <= 70
                    global Bd = params.boundary_layer_height_3(Ht[pos2idx(i,j)],T[pos2idx(i,j)]) #boundary layer, m
                else 
                    global Bd = params.boundary_layer_height_4 #boundary layer, m
                end

                #salinity
                Sal[pos2idx(i,j)] = S(T[pos2idx(i,j)],i,Vavg(Ht[pos2idx(i,j)],Bd)) #kg/m3
                #height 
                Ht[pos2idx(i,j)] = Hght(i,T[pos2idx(i,0)],Sal[pos2idx(i,j)]) #m
                #increments in z direction
                dz_v[pos2idx(i,j)] = dz(i,T[pos2idx(i,0)],Sal[pos2idx(i,j)]) #m
            
                if j >= 4*Nz/5 #bottom 20% of channel
                    V_profile[pos2idx(i,j)] = params.velocity_profile_nw(dz_v[pos2idx(i,j)]*j, Ht[pos2idx(i,j)], Bd) #non-wake flow, m/hr
                else
                    V_profile[pos2idx(i,j)] = params.velocity_profile_w(dz_v[pos2idx(i,j)]*j,Ht[pos2idx(i,j)], Bd) #wake flow in rest of channel, m/hr
                end
                

            else
                Sal[pos2idx(i,j)] = S(T[pos2idx(i,j)],i,Vavg_lam(Ht[pos2idx(i,j)],T[pos2idx(i,j)]))
               
                dz_v[pos2idx(i,j)] = dz(i,T[pos2idx(i,0)],Sal[pos2idx(i,j)])
                Ht[pos2idx(i,j)] = Hght(i,T[pos2idx(i,0)],Sal[pos2idx(i,j)])
                   
                V_profile[pos2idx(i,j)] =  params.velocity_profile_lam(dz_v[pos2idx(i,j)]*j,Ht[pos2idx(i,j)],T[pos2idx(i,j)])
                Bd = 0
                 
            end
        end
    end
    


    #This function defines the heat flux from direct solar radiation

    Q_Solar = (1-0.03) * alpha_w * GHI * 3600;                                  #goes in, positive. removed * (1-f_a) as the fraction of biomass converting light into chemical energy is very small

    # Source: Yadala and Cremaschi, 2016
    
    #This function defines the heat flux from reradiation from the pond (Q_out)

    Q_Rerad(Temp) = -1.0 * sigma * epsilon_water .* Temp.^4 * 3600              #goes out, negative

    #This function defines the heat flux from longwave atmospheric radiation
    Q_Longwave_Atmo = epsilon_water * epsilon_air * sigma .* Tamb.^4 * 3600     #goes in, positive

    #This function calculates the evaporation rate of the pond as well as the
    #cooling effect due to that evaporation. Straight from Yadala and Cremaschi, 2016.

    #Calculate the Reynold's Number
    kin_visc_a = dyn_visc_air / rho_air      # m^2/s
    L_c = L #m, characteristic length
    Re_L = L_c * WNDSPD / kin_visc_a #Reynold's number (air over water)

    Q_Evap(Temp) = params.evaporation_heat_flux(Temp,WNDSPD,RH,P_a)


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

    ##T(-1,j) = T(0,j)
    ##T(Ny+1,j) = T(Ny,j)

    
  

    for i = 1:Ny-1
       
 
        # BC 2
        dT[pos2idx(i, 0)] = (  k_water*(T[pos2idx(i-1, 0)] - 2*T[pos2idx(i, 0)] + T[pos2idx(min(Ny,i+1), 0)]) / dy^2    #conduction in y-dir
                             + k_water*(T[pos2idx(i, 0)] + T[pos2idx(i,0+1)] - 2*T[pos2idx(i, 0)]) / dz_v[pos2idx(i,0)]^2               #conduction in z-dir at boundary
                             + Q_sum1(T[pos2idx(i,0)]) * dy * W                                                         #heat generation at boundary
                             - V_profile[pos2idx(i,0)] * ( T[pos2idx(i,0)] - T[pos2idx(i-1,0)]) / dy                                       #heat convection at boundary
                             ) / (rho_water * cp_water)

        # BC 3
        dT[pos2idx(i, Nz)] = (  k_water*(T[pos2idx(i-1, Nz)] - 2*T[pos2idx(i, Nz)] + T[pos2idx(min(Ny,i+1), Nz)]) / dy^2 #conduction in y-dir
                              + k_water*(T[pos2idx(i, Nz-1)] + T[pos2idx(i,Nz)] - 2*T[pos2idx(i, Nz)]) / dz_v[pos2idx(i,Nz)]^2            #conduction in z-dir at boundary
                              + Q_sum2(T[pos2idx(i,Nz)]) * dy * W                                                        #heat generation at boundary
                              - V_profile[pos2idx(i,Nz)]*(T[pos2idx(i,Nz)] - T[pos2idx(i-1,Nz)])/ dy                                        #heat convection at boundary
                             ) / (rho_water*cp_water)
    end

    for j = 0:Nz

        dT[pos2idx(0, j)] = 0 

        dT[pos2idx(Ny, j)] = (  k_water * (T[pos2idx(Ny, max(0,j-1))] - 2*T[pos2idx(Ny,j)] + T[pos2idx(Ny, min(Nz,j+1))]) / dz_v[pos2idx(Ny,j)]^2    #conduction in z-dir
                              + k_water * (T[pos2idx(Ny-1, j)] - 2*T[pos2idx(Ny,j)] + T[pos2idx(Ny,j)]) / dy^2                      #conduction in y-dir at boundary
                              - V_profile[pos2idx(Ny,j)] / dy                                      #heat convection at boundary
                             ) / (rho_water*cp_water)
    end

    for i=1:Ny-1
        for j=1:Nz-1
            dT[pos2idx(i,j)] = (   k_water * (T[pos2idx(i-1, j)] - 2*T[pos2idx(i, j)] + T[pos2idx(i+1, j)]) / dy^2      #conduction in y-dir
                                 + k_water * (T[pos2idx(i, j-1)] - 2*T[pos2idx(i, j)] + T[pos2idx(i, j+1)]) / dz_v[pos2idx(i,j)]^2      #conduction in z-dir
                                 - V_profile[pos2idx(i,j)] / dy                         #heat convection
                               ) / (rho_water*cp_water)
        end
    end


    #@show dT
    @views dX[Nelements+1:2*Nelements] .= dT        # order matters!  The @views operator takes a slice out of an array without making a copy.
    nothing
end
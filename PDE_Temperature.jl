
function HeatTransfer!(dX, T, params, t)
    #Units for params can be found in LoadParameters.jl

    #Environmental Data
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data
    data_begin = params.data_begin

    #Timestep
    t_hour1 = floor(Int64, t)
    t_hour2 = floor(Int64, t)+1

    #Adjusts Environmental Data for each timestep
    GHI = GHI_Data[data_begin + t_hour1] * (t-t_hour1) + GHI_Data[data_begin + t_hour2] * (t_hour2 - t)
    WNDSPD = max.((WNDSPD_Data[data_begin + t_hour1] * (t-t_hour1) + WNDSPD_Data[data_begin + t_hour2] * (t_hour2 - t)),0)
    RH = RH_Data[data_begin + t_hour1] * (t-t_hour1) + RH_Data[data_begin + t_hour2] * (t_hour2 - t)
    Tamb = max.(((Tamb_Data[data_begin + t_hour1] * (t-t_hour1) + Tamb_Data[data_begin + t_hour2] * (t_hour2 - t))),0)

    #Water and Air Properties
    rho_water = params.density_water(Tamb)
    rho_air = params.density_air(Tamb)
    dyn_visc_air = params.dynamic_viscosity_air(Tamb)
    cp_water = params.specific_heat_capacity_water(Tamb)
    pr_air = params.prandtl_air(Tamb)
    P_a = params.saturated_vapor_pressure_water_air(Tamb)
    k_water = params.thermal_conductivity_water(Tamb)
    k_air = params.thermal_conductivity_air(Tamb)
    epsilon_water = params.emissivity_water
    epsilon_air = params.emissivity_air
    alpha_w = params.solar_reflectance_water

    #Ground Properties
    diff_concrete = params.thermal_diffusivity_concrete
    ground_temperature = params.ground_temperature
    k_concrete = params.thermal_conductivity_concrete

     #Geometric Properties
     L = params.reactor_length
     W = params.reactor_width
     H = params.reactor_initial_liquid_level
     Ny = params.num_odes_y
     Nz = params.num_odes_z
     dy = L / Ny

     #height and height increment as variables of distance across length of reactor
     Hght(x,S) = params.height(Tamb,S,RH,WNDSPD,P_a,x)
     dz(x,S) = Hght(x,S)/Nz
 
     #Reindexing of 2D position in reactor to be contained in 1D vector
     pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
     idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
 
     #Number of total elements in vector
     Nelements = (Ny+1) * (Nz+1) #number of elements in vectors containing info in y and z-direction
     Nelements1 = Ny + 1 #number of elements in vectors containing info in y-direction
   
     T = max.(T, 0.0) #don't allow negative values
     dT = zeros( Nelements, 1) #dT/dt
 
     #Other parameters
     P_a = params.saturated_vapor_pressure_water_air(Tamb) 
     M_Evap = params.evaporation_mass_flux(Tamb,WNDSPD,RH,P_a)
     sigma = params.stefan_boltzmann_constant
     S(x) = params.salinity(Tamb,RH,WNDSPD,P_a,x)
     Re(H,T,V) = params.reynolds_number(H,T,V)
 
     #Vavg determined assuming conservation of kinetic energy (evaporated water has 0 kinetic energy)
     M = zeros(Nelements1) #mass of water in vertical slice (kg)
     Vavg = zeros(Nelements1) #average velocity

     #Initialize vectors for vertical increments, salinity, height
     dz_v = zeros(Nelements1)
     Sal = zeros(Nelements1)
     Ht = zeros(Nelements1)
 
     #Initialize Mass and Avg Velocity
     M[pos2idx(0,0)] = W*dy*H*rho_water
     Vavg[pos2idx(0,0)] = params.volumetric_flow_rate_o/(W*H)
 
     for i in 1:Ny
         M[pos2idx(i,0)] = M[pos2idx(i-1,0)] - dy*(1/Vavg[pos2idx(i-1,0)])*M_Evap*dy*W
         Vavg[pos2idx(i,0)] = sqrt((M[pos2idx(i-1,0)]*(Vavg[pos2idx(i-1,0)])^2)/(M[pos2idx(i,0)]))
     end

    for i in 0:Ny
        #salinity
        Sal[pos2idx(i,0)] = S(i) #kg/m3
        #height 
        Ht[pos2idx(i,0)] = Hght(i,Vavg[pos2idx(i,0)]) #m
        #increments in z direction
        dz_v[pos2idx(i,0)] = dz(i,Vavg[pos2idx(i,0)]) #m
    end
 
     V_profile = zeros(Nelements, 1) #velocity profile (velocity at different points across vertical axis)
 
     for i in 0:Ny
         for j in 0:Nz
 
                 #different velocity profiles under laminar and turbulent conditions
                 if Re(Ht[pos2idx(0,0)],T[pos2idx(0,0)],Vavg[pos2idx(0,0)]) > 2000 #turbulent
                     V_profile[pos2idx(i,j)] = Vavg[pos2idx(i,0)]
                 else #laminar
                     V_profile[pos2idx(i,j)] = params.velocity_profile_lam(Vavg[pos2idx(i,0)],j,Ht[pos2idx(i,0)])
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

    Q_Evap(Temp) = -params.evaporation_heat_flux(Temp,WNDSPD,RH,P_a) #goes out, negative


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
   
    #BC1: T(y,z) at y = 0 is constant [Tin]
    #BC2: dT(y,z) at z = 0 includes the d2T/dy2, Vy*dT/dy, and Q terms (Qsol, Qrerad, Qlw, Qcov, Qevap)
    #BC3: dT(y,z) at z = H includes the d2T/dy2, Vy*dT/dy, and Q terms (Qground)

    ##T(-1,j) = T(0,j)
    ##T(Ny+1,j) = T(Ny,j)


    for i = 1:Ny
       
 
        # BC 2
        dT[pos2idx(i, 0)] = (  (k_water*(T[pos2idx(i-1, 0)] - 2*T[pos2idx(i, 0)] + T[pos2idx(min(Ny,i+1), 0)]) / dy^2    #conduction in y-dir
                             + Q_sum1(T[pos2idx(i,0)]) * dy * W)                                                        #heat generation at boundary
                             - (V_profile[pos2idx(i,0)]*(T[pos2idx(i,0)] - T[pos2idx(i-1,0)])/dy)                                       #heat convection at boundary
                             ) / (rho_water * cp_water)

        # BC 3
        dT[pos2idx(i, Nz)] = (  (k_water*(T[pos2idx(i-1, Nz)] - 2*T[pos2idx(i, Nz)] + T[pos2idx(min(Ny,i+1), Nz)]) / dy^2 #conduction in y-dir
                              + Q_sum2(T[pos2idx(i,Nz)]) * dy * W)                                                       #heat generation at boundary
                              - (V_profile[pos2idx(i,Nz)]*(T[pos2idx(i,Nz)] - T[pos2idx(i-1,Nz)])/dy)                                       #heat convection at boundary
                             ) / (rho_water*cp_water)
    end

    for j = 0:Nz

        dT[pos2idx(0, j)] = 0 

    end

    for i=1:Ny
        for j=1:Nz-1
            dT[pos2idx(i,j)] = (   (k_water * (T[pos2idx(i-1, j)] - 2*T[pos2idx(i, j)] + T[pos2idx(min(Ny,i+1), j)]) / dy^2      #conduction in y-dir
                                 + k_water * (T[pos2idx(i, j-1)] - 2*T[pos2idx(i, j)] + T[pos2idx(i, j+1)]) / dz_v[pos2idx(i,0)]^2)     #conduction in z-dir
                                 - (V_profile[pos2idx(i,j)]*(T[pos2idx(i,j)] - T[pos2idx(i-1,j)])/dy)                         #heat convection
                               ) / (rho_water*cp_water)
        end
    end

    @views dX[Nelements+1:2*Nelements] .= dT        # order matters!  The @views operator takes a slice out of an array without making a copy.
    nothing
end

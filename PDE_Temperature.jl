
function HeatTransfer!(dX, T, A, params, t)
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
    rho_water(T) = params.density_water(T)
    rho_air = params.density_air(Tamb)
    dxn_visc_air = params.dynamic_viscosity_air(Tamb)
    cp_water(T) = params.specific_heat_capacity_water(T)
    pr_air = params.prandtl_air(Tamb)
    P_a = params.saturated_vapor_pressure_water_air(Tamb)
    k_air = params.thermal_conductivity_air(Tamb)
    k_water(T) = params.thermal_conductivity_water(T)
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
     Nx = params.num_odes_x
     Ny = params.num_odes_y
     dx = L / Nx

     #height and height increment as variables of distance across length of reactor
     Hght(x,S) = params.height(Tamb,S,RH,WNDSPD,P_a,x)
     dy(x,S) = Hght(x,S)/Ny
 
     #Reindexing of 2D position in reactor to be contained in 1D vector
     pos2idx(y) = (y.+1)
     idx2pos(pos) = [Integer(pos - 1 - (Nx+1) * floor( (pos-1) ./ (Nx.+1))), Integer(floor( (pos-1) ./ (Nx.+1)))]
 
     #Number of total elements in vector
     Nelements = (Nx+1) #number of elements in vectors containing info in y and z-direction
   
     T = max.(T, 0.0) #don't allow negative values
     dT = zeros( Nelements, 1) #dT/dt
     dA = zeros(Nelements)
 
     #Other parameters
     P_a = params.saturated_vapor_pressure_water_air(Tamb) 
     sigma = params.stefan_boltzmann_constant
     Re(H,T,V) = params.reynolds_number(H,T,V)
     alpha = params.alpha
     m = params.m
     Vavg = zeros(Nelements) #average velocity

     #Initialize vectors for vertical increments, salinity, height
     H = zeros(Nelements)

    for i in 0:Nx
        Vavg[pos2idx(i)] = alpha*A[pos2idx(i)]^(m-1)
        H[pos2idx(i)] = A[pos2idx(i)]/(W)
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
    kin_visc_a = dxn_visc_air / rho_air      # m^2/s
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
    #BC2: dT(y,z) at z = 0 includes the d2T/dx2, Vy*dT/dx, and Q terms (Qsol, Qrerad, Qlw, Qcov, Qevap)
    #BC3: dT(y,z) at z = H includes the d2T/dx2, Vy*dT/dx, and Q terms (Qground)

    ##T(-1,j) = T(0,j)
    ##T(Nx+1,j) = T(Nx,j)


    T_conv = zeros(Nx+1)
    T_str = zeros(Nx+1)
    strip_Q = params.volumetric_flow_rate_strip
    Term_adj = zeros(Nx+1)
    Terms_flux = zeros(Nx+1)
    Terms_conduction = zeros(Nx+1)

    dT[pos2idx(0)] = 0 

    q(Temp) = params.lat_flow(Temp,WNDSPD,RH,P_a)



    for i=0:Nx
        dA[pos2idx(i)] = q(T[pos2idx(i)]) -alpha*m*(((A[pos2idx(i)] + A[pos2idx(max(i-1,0))])/2)^(m-1))*(A[pos2idx(i)]-A[pos2idx(max(i-1,0))])/dx
    end

    for i=1:Nx
            T_conv[pos2idx(i)] = -(A[pos2idx(i)]*Vavg[pos2idx(i)]*T[pos2idx(i)] - A[pos2idx(i-1)]*Vavg[pos2idx(i-1)]*T[pos2idx(i-1)])/(A[pos2idx(i)]*dx) #K/hr
            T_str[pos2idx(i)] = (strip_Q*params.input_temperature)/(A[pos2idx(i)]*L) #K/hr
            Term_adj[pos2idx(i)] = -dA[pos2idx(i)]*T[pos2idx(i)]/(A[pos2idx(i)]) #K/hr
            Terms_flux[pos2idx(i)] = (Q_sum1(T[pos2idx(i)])+Q_sum2(T[pos2idx(i)]))/(rho_water(T[pos2idx(i)])*cp_water(T[pos2idx(i)])*(H[pos2idx(i)]))
            Terms_conduction[pos2idx(i)] = (k_water(T[pos2idx(i)]) * (T[pos2idx(i-1)] - 2*T[pos2idx(i)] + T[pos2idx(min(Nx,i+1))]) / dx^2)/(rho_water(T[pos2idx(i)])*cp_water(T[pos2idx(i)]))

            dT[pos2idx(i)] = (( (k_water(T[pos2idx(i)]) * (T[pos2idx(i-1)] - 2*T[pos2idx(i)] + T[pos2idx(min(Nx,i+1))]) / dx^2) 
                               ) / (rho_water(T[pos2idx(i)])*cp_water(T[pos2idx(i)]))#K/hr
                               + Terms_flux[pos2idx(i)] #K/hr
                               + T_str[pos2idx(i)] 
                               + T_conv[pos2idx(i)]
                               + Term_adj[pos2idx(i)])
            
                               
    end

    @views dX[Nelements+1:2*Nelements] .= dT        # order matters!  The @views operator takes a slice out of an array without making a copy.
    nothing
end

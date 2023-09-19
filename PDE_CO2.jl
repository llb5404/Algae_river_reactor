
function PDE_CO2!(dX, C, C_biomass, Temperature, params, t)
    #Units for params can be found in LoadParameters.jl
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data

    t_hour1 = floor(Int64, t)
    t_hour2 = floor(Int64, t)+1

    data_begin = params.data_begin

    GHI = GHI_Data[data_begin + t_hour1] * (t-t_hour1) + GHI_Data[data_begin + t_hour2] * (t_hour2 - t)
    Tamb = max.(((Tamb_Data[data_begin + t_hour1] * (t-t_hour1) + Tamb_Data[data_begin + t_hour2] * (t_hour2 - t))),0)
    RH = RH_Data[data_begin + t_hour1] * (t-t_hour1) + RH_Data[data_begin + t_hour2] * (t_hour2 - t)
    WNDSPD = max.((WNDSPD_Data[data_begin + t_hour1] * (t-t_hour1) + WNDSPD_Data[data_begin + t_hour2] * (t_hour2 - t)),0)

    mu_model(T,S,M_biomass_z,z,dz,C_co2) = params.biomass_specific_growth_rate(T, S, M_biomass_z, GHI, z, dz, C_co2)
    co2_per_biomass = params.co2_per_biomass
    L = params.reactor_length
    W = params.reactor_width
    H = params.reactor_initial_liquid_level
    Ny = params.num_odes_y
    Nz = params.num_odes_z

    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    dy = L / Ny

    D_co2 = params.diffusion_coeff_co2_water    #m^2 / sec, depends on Temperature
    S_co2 = params.solubility_co2_water         #mole fraction of dissolved CO2 in water at equilibrium, depends on Temperature.
    density_water = params.density_water        # kg / m^3, depends on Temperature.
    molecular_weight_co2 = params.molecular_weight_co2
    molecular_weight_water = params.molecular_weight_water
    P_a = params.saturated_vapor_pressure_water_air(Tamb) 

    Nelements = (Ny+1) * (Nz+1)
    C = max.(C, 0.0) #don't allow negative values
    dCO2 = zeros( Nelements, 1)
    H_o = params.reactor_initial_liquid_level
    K(T,S) = params.dVavgdx(T,S,RH,WNDSPD,P_a)
    M_Evap(T) = params.evaporation_mass_flux(T,WNDSPD,RH,P_a)
    rho_solution(S) = params.density_solution(Tamb,S)
    S(T,x,V) = params.salinity(T,WNDSPD,RH,P_a,x,V)

   

    Hght(x,T,S) = params.height(x,T,S,WNDSPD,RH,P_a,H_o)
    dz(x,T,S) = Hght(x,T,S)/Nz

    ##Velocity Profile
    Re(H,T,V) = params.reynolds_number(H,T,V)
    Re_star(H,T) = params.rough_reynolds_number(H,T)
    V_profile = zeros(Nelements, 1)
    Vavg_lam(H,T) = params.average_flow_velocity_lam(H,T)
    Vavg(H,B) = params.average_flow_velocity(H,B)
    dz_v = zeros(Nelements,1)
    Sal = zeros(Nelements, 1)
    Ht = zeros(Nelements,1)
    R_co2 = zeros(Nelements,1)
    mu_v = zeros(Nelements, 1)


    for i in 0:Ny
        for j in 0:Nz
            #cut-off for laminar flow = 500
            if Re(Hght(i,Temperature[pos2idx(i,j)],S(Temperature[pos2idx(i,j)],i,V_profile[pos2idx(i,j)])),Temperature[pos2idx(i,j)],Vavg_lam(Ht[pos2idx(i,j)],Temperature[pos2idx(i,j)])) > 500
                #Turbulent Flow: choose boundary layer thickness based on "rough reynolds number"
                if Re_star(Hght(i,Temperature[pos2idx(i,j)],S(Temperature[pos2idx(i,j)],i,V_profile[pos2idx(i,j)])),Temperature[pos2idx(i,j)]) <= 4
                    global Bd = params.boundary_layer_height_1(Ht[pos2idx(i,j)],Temperature[pos2idx(i,j)]) #boundary layer, m
                elseif Re_star(Hght(i,Temperature[pos2idx(i,j)],S(Temperature[pos2idx(i,j)],i,V_profile[pos2idx(i,j)])),Temperature[pos2idx(i,j)]) <= 11
                    global Bd = params.boundary_layer_height_2(Ht[pos2idx(i,j)],Temperature[pos2idx(i,j)]) #boundary layer, m
                elseif Re_star(Hght(i,Temperature[pos2idx(i,j)],S(Temperature[pos2idx(i,j)],i,V_profile[pos2idx(i,j)])),Temperature[pos2idx(i,j)]) <= 70
                    global Bd = params.boundary_layer_height_3(Ht[pos2idx(i,j)],Temperature[pos2idx(i,j)]) #boundary layer, m
                else 
                    global Bd = params.boundary_layer_height_4 #boundary layer, m
                end

                #salinity
                Sal[pos2idx(i,j)] = S(Temperature[pos2idx(i,j)],i,Vavg(Ht[pos2idx(i,j)],Bd)) #kg/m3
                #height 
                Ht[pos2idx(i,j)] = Hght(i,Temperature[pos2idx(i,0)],Sal[pos2idx(i,j)]) #m
                #increments in z direction
                dz_v[pos2idx(i,j)] = dz(i,Temperature[pos2idx(i,0)],Sal[pos2idx(i,j)]) #m
            
                if j >= 4*Nz/5 #bottom 20% of channel
                    V_profile[pos2idx(i,j)] = params.velocity_profile_nw(dz_v[pos2idx(i,j)]*j, Ht[pos2idx(i,j)], Bd) #non-wake flow, m/hr
                else
                    V_profile[pos2idx(i,j)] = params.velocity_profile_w(dz_v[pos2idx(i,j)]*j,Ht[pos2idx(i,j)], Bd) #wake flow in rest of channel, m/hr
                end
                

            else
                Sal[pos2idx(i,j)] = S(Temperature[pos2idx(i,j)],i,Vavg_lam(Ht[pos2idx(i,j)],Temperature[pos2idx(i,j)]))
               
                dz_v[pos2idx(i,j)] = dz(i,Temperature[pos2idx(i,0)],Sal[pos2idx(i,j)])
                Ht[pos2idx(i,j)] = Hght(i,Temperature[pos2idx(i,0)],Sal[pos2idx(i,j)])
                   
                V_profile[pos2idx(i,j)] =  params.velocity_profile_lam(dz_v[pos2idx(i,j)]*j,Ht[pos2idx(i,j)],Temperature[pos2idx(i,j)])
                Bd = 0
                 
            end
        end
    end

    #Vector of specific growth rate values
    mu_v = zeros(Nelements, 1)
    for i in 0:Ny
        for j in 0:Nz
            mu_v[pos2idx(i,j)] = mu_model(Temperature[pos2idx(i,j)],Sal[pos2idx(i,j)],C_biomass[pos2idx(i,0:Nz)],j,dz_v[pos2idx(i,j)],C[pos2idx(i,j)]) #1/hr
        end
    end

    # CO2 per biomass x biomass production rate [kg CO2 / kg biomass  x  kg biomass/hour = kg CO2 / hour]
    R_co2 = zeros(Nelements,1)
    for i in 1:Ny
      for j in 1:Nz
        R_co2[pos2idx(i,j)] = co2_per_biomass * mu_v[pos2idx(i,j)] * C_biomass[j]/(dz_v[pos2idx(i,j)]*dy*W)
      end
    end
    
    

    # ==========  DISSOLVED CARBON DIOXIDE BALANCE ========
    # C is dissolved CO2 [kg CO2 / m^3 water ]

    #BC1: CO2(y,z)  at y = 0 is constant [CO2_in]
    #     dCO2(y,z) at y = 0 is zero.

    for i = 1:Ny
      C[pos2idx(i,0)] =  (S_co2(Temperature[pos2idx(i,0)])* density_water(Temperature[pos2idx(i,0)])* molecular_weight_co2/ molecular_weight_water)
    end

    dCO2[pos2idx(1:Ny,0)] .= 0.0

    #y = dy .. Ny,  z = 0 or z = Nz
    for i=1:Ny-1

        dCO2[pos2idx(i,Nz)] =  ( + D_co2(Temperature[pos2idx(i,Nz)]) * (C[pos2idx(i-1,Nz)] + C[pos2idx(i+1,Nz)] - 2*C[pos2idx(i,Nz)]) / dy^2 #small     #diffusion CO2 in y-direction
                                 + D_co2(Temperature[pos2idx(i,Nz)]) * (C[pos2idx(i,Nz-1)] + C[pos2idx(i,Nz)] - 2*C[pos2idx(i,Nz)]) / dz_v[pos2idx(i,Nz)]^2              #diffusion CO2 in z-direction
                                 - V_profile[pos2idx(i,Nz)] * (C[pos2idx(i,Nz)] - C[pos2idx(i-1,Nz)]) / dy                                         #convection of CO2 at bottom (zero)
                                 - R_co2[pos2idx(i,Nz)]                                                    #consumption of CO2 by biomass growth
                               )
    end

    dCO2[pos2idx(Ny,Nz)] =     ( + D_co2(Temperature[pos2idx(Ny,Nz)]) * (C[pos2idx(Ny-1,Nz)] + C[pos2idx(Ny,Nz)] - 2*C[pos2idx(Ny,Nz)]) / dy^2 #small       #diffusion CO2 in y-direction
                                 + D_co2(Temperature[pos2idx(Ny,Nz)]) * (C[pos2idx(Ny,Nz-1)] + C[pos2idx(Ny,Nz)] - 2*C[pos2idx(Ny,Nz)]) / dz_v[pos2idx(Ny,Nz)]^2              #diffusion CO2 in z-direction
                                 - V_profile[pos2idx(Ny,Nz)] * (C[pos2idx(Ny,Nz)] - C[pos2idx(Ny-1,Nz)]) / dy                                         #convection of CO2 at bottom (zero)
                                 - R_co2[pos2idx(Ny,Nz)]                                                  #consumption of CO2 by biomass growth
                               )

    for i=1:Ny-1
        for j=1:Nz-1
            dCO2[pos2idx(i,j)]=  ( + D_co2(Temperature[pos2idx(i,j)]) * (C[pos2idx(i-1,j)] + C[pos2idx(i+1,j)] - 2*C[pos2idx(i,j)]) / dy^2             #diffusion CO2 in y-direction
                                   + D_co2(Temperature[pos2idx(i,j)]) * (C[pos2idx(i,j-1)] + C[pos2idx(i,j+1)] - 2*C[pos2idx(i,j)]) / dz_v[pos2idx(i,j)]^2             #diffusion CO2 in z-direction
                                   - V_profile[pos2idx(i,j)] * (C[pos2idx(i,j)] - C[pos2idx(i-1,j)]) / dy                           #convection of CO2 in y-direction
                                   - R_co2[pos2idx(i,j)]                                                 #consumption of CO2 by biomass growth
                                 )

        end
    end

    @views dX[1+2*Nelements:3*Nelements] .= dCO2        # order matters! The @views operator takes a slice out of an array without making a copy.
    nothing
end
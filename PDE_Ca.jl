function PDE_Ca!(dX, C, CO3, Temperature, params, t)
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
 
    molecular_weight_water = params.molecular_weight_water
    P_a = params.saturated_vapor_pressure_water_air(Tamb)
    
    #density to M conversions
    molecular_weight_co2 = params.molecular_weight_co2
    co2_to_M = (1/(molecular_weight_co2*1000)) #converts from kg/m3 to mol/L: mol*m3/kg-L
    molecular_weight_h2co3 = params.molecular_weight_h2co3
    h2co3_to_M = (1/(molecular_weight_h2co3*1000))
    molecular_weight_hco3 = params.molecular_weight_hco3
    hco3_to_M = (1/(molecular_weight_hco3*1000))
    molecular_weight_co3 = params.molecular_weight_co3
    co3_to_M = (1/(molecular_weight_co3*1000))
    molecular_weight_water = params.molecular_weight_water
    molecular_weight_ca = params.molecular_weight_ca
    ca_to_M = (1/(molecular_weight_ca*1000))
    molecular_weight_caco3 = params.molecular_weight_caco3
    caco3_to_M = (1/(molecular_weight_caco3*1000))

    #initial conditions
    pH_o = params.initial_pH
    HCO3_o = params.initial_hco3
    CO3_o = params.initial_co3
    Ca_o = params.initial_ca
    CaCO3_o = params.initial_caco3

    #kinetics
    k1 = params.k_1
    kneg1 = params.k_neg1
    k2 = params.k_2
    kneg2 = params.k_neg2
    k3 = params.k_3
    kneg3 = params.k_neg3
    k4 = params.k_4
    kneg4 = params.k_neg4
    k5 = params.k_5
    kneg5 = params.k_neg5
    sub = params.substrate

    Nelements = (Ny+1) * (Nz+1)
    C = max.(C, 0.0) #don't allow negative values
    dCa = zeros( Nelements, 1)
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

    Z = zeros(Nelements,1)
    for i in 0:Ny
        for j in 0:Nz
            Z[pos2idx(i,j)] = CaCO3_o*caco3_to_M - C[pos2idx(i,j)]*ca_to_M + Ca_o*ca_to_M
        end
    end

   
    # ==========  Ca BALANCE ========


    for i=1:Ny

        dCa[pos2idx(i,Nz)] =  (  (-k5*CO3[pos2idx(i,Nz)]*C[pos2idx(i,Nz)]*co3_to_M*ca_to_M + kneg5*min(sub,Z[pos2idx(i,Nz)]))*(1/ca_to_M) #reaction kinetics of H2CO3
                                 - V_profile[pos2idx(i,Nz)] * (C[pos2idx(i,Nz)] - C[pos2idx(i-1,Nz)]) / dy    #convection of H2CO3
                               )
    end


    for i=1:Ny
        for j=0:Nz-1
            dCa[pos2idx(i,Nz)] =  (  (-k5*CO3[pos2idx(i,j)]*C[pos2idx(i,j)]*co3_to_M*ca_to_M + kneg5*min(sub,Z[pos2idx(i,j)]))*(1/ca_to_M) #reaction kinetics of H2CO3
                                 - V_profile[pos2idx(i,j)] * (C[pos2idx(i,j)] - C[pos2idx(i-1,j)]) / dy    #convection of H2CO3
                               )

        end
    end

    
    
    @views dX[1+6*Nelements:7*Nelements] .= dCa        # order matters! The @views operator takes a slice out of an array without making a copy.
    nothing
end
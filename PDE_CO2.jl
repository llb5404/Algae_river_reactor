
function PDE_CO2!(dX, C, Height, C_biomass, Temperature, params, t)

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

    Isat = params.max_biomass_light_saturation
    Cmax = params.max_biomass_concentration
    mu_model = params.biomass_specific_growth_rate
    co2_per_biomass = params.co2_per_biomass
    Q = params.volumetric_flow_rate
    V_profile(z,H) = params.velocity_profile(z,H,Tamb)
    L = params.reactor_length
    W = params.reactor_width
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

    Hght(x,T,S) = params.height(x,T,S,WNDSPD,RH,P_a,H_o)
    dz(x,T,S) = Hght(x,T,S)/Nz
    

    mu(Mz, z, dz_y, C_co2) = mu_model(Mz, GHI, z, dz_y, C_co2) # biomass specific growth rate, 1 / hour

    R_co2(Mz, z, dz_y, C_co2) = co2_per_biomass * mu(Mz, z, dz_y, C_co2)* Mz[z]    # CO2 per biomass x biomass production rate [kg CO2 / kg biomass  x  kg biomass/hour = kg CO2 / hour]
  
    # ==========  DISSOLVED CARBON DIOXIDE BALANCE ========
    # C is dissolved CO2 [kg CO2 / m^3 water ]

    #BC1: CO2(y,z)  at y = 0 is constant [CO2_in]
    #     dCO2(y,z) at y = 0 is zero.
    dCO2[pos2idx(0,0:Nz)] .= 0.0

    #BC2: Constant CO2(y,z) at z = 0 (surface).  Assume equilibrium at water-air interface. CO2 concentration is determined by solubility.
    # S_co2 [mol CO2 / mol water] * density_water [kg water / m^3 water] * mw_co2 [kg CO2 / mole CO2] / mw_water [kg water / mole water] => kg CO2 / m^3 water
    
    for i = 1:Ny
      C[pos2idx(i,0)] =  (S_co2(Temperature[pos2idx(i,0)])* density_water(Temperature[pos2idx(i,0)])* molecular_weight_co2/ molecular_weight_water)
    end

    dCO2[pos2idx(1:Ny,0)] .= 0.0

    #y = dy .. Ny,  z = 0 or z = Nz
    for i=1:Ny-1

        dCO2[pos2idx(i,Nz)] =  ( + D_co2(Temperature[pos2idx(i,Nz)]) * (C[pos2idx(i-1,Nz)] + C[pos2idx(i+1,Nz)] - 2*C[pos2idx(i,Nz)]) / dy^2 #small     #diffusion CO2 in y-direction
                                 + D_co2(Temperature[pos2idx(i,Nz)]) * (C[pos2idx(i,Nz-1)] + C[pos2idx(i,Nz)] - 2*C[pos2idx(i,Nz)]) / dz(i,Temperature[pos2idx(i,Nz)],0)^2              #diffusion CO2 in z-direction
                                 - V_profile(dz(i,Temperature[pos2idx(i,Nz)],0)*Nz, Hght(i,Temperature[pos2idx(i,Nz)],0)) * (C[pos2idx(i,Nz)]/ - C[pos2idx(i-1,Nz)]) / dy                                         #convection of CO2 at bottom (zero)
                                 - R_co2(M_biomass[pos2idx(i,0:Nz)], Nz, dz(i,Temperature[pos2idx(i,Nz)],0), C[pos2idx(i,Nz)])                                                    #consumption of CO2 by biomass growth
                               )
    end

    dCO2[pos2idx(Ny,Nz)] =     ( + D_co2(Temperature[pos2idx(Ny,Nz)]) * (C[pos2idx(Ny-1,Nz)] + C[pos2idx(Ny,Nz)]/ - 2*C[pos2idx(Ny,Nz)]) / dy^2 #small       #diffusion CO2 in y-direction
                                 + D_co2(Temperature[pos2idx(Ny,Nz)]) * (C[pos2idx(Ny,Nz-1)] + C[pos2idx(Ny,Nz)] - 2*C[pos2idx(Ny,Nz)]) / dz(Ny,Temperature[pos2idx(Ny,Nz)],0)^2              #diffusion CO2 in z-direction
                                 - V_profile(dz(Ny,Temperature[pos2idx(Ny,Nz)],0)*Nz, Hght(Ny,Temperature[pos2idx(Ny,Nz)],0)) * (C[pos2idx(Ny,Nz)]/ - C[pos2idx(Ny-1,Nz)]) / dy                                         #convection of CO2 at bottom (zero)
                                 - R_co2(M_biomass[pos2idx(Ny,0:Nz)], Nz, dz(Ny,Temperature[pos2idx(Ny,Nz)],0), C[pos2idx(Ny,Nz)])                                                    #consumption of CO2 by biomass growth
                               )

    for i=1:Ny-1
        for j=1:Nz-1
            dCO2[pos2idx(i,j)]=  ( + D_co2(Temperature[pos2idx(i,j)]) * (C[pos2idx(i-1,j)] + C[pos2idx(i+1,j)] - 2*C[pos2idx(i,j)]) / dy^2             #diffusion CO2 in y-direction
                                   + D_co2(Temperature[pos2idx(i,j)]) * (C[pos2idx(i,j-1)] + C[pos2idx(i,j+1)] - 2*C[pos2idx(i,j)]) / dz(i,Temperature[pos2idx(i,j)],0)^2             #diffusion CO2 in z-direction
                                   - V_profile(dz(i,Temperature[pos2idx(i,j)],0)*j, Hght(i,Temperature[pos2idx(i,j)],0)) * (C[pos2idx(i,j)]/ - C[pos2idx(i-1,j)]) / dy                           #convection of CO2 in y-direction
                                   - R_co2(M_biomass[pos2idx(i,0:Nz)], j, dz(i,Temperature[pos2idx(i,j)],0), C[pos2idx(i,j)])                                                   #consumption of CO2 by biomass growth
                                 )
        end
    end

    @views dX[1+2*Nelements:3*Nelements] .= dCO2        # order matters! The @views operator takes a slice out of an array without making a copy.
    nothing
end
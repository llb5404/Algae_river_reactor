
function PDE_CO2!(dX, C, DIC, C_biomass, Temperature, params, t)
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

    mu_model(T,S,Co2) = params.biomass_specific_growth_rate(T, S, Co2)
    co2_per_biomass = params.co2_per_biomass
    bm = params.max_biomass_specific_growth_rate
    L = params.reactor_length
    W = params.reactor_width
    H = params.reactor_initial_liquid_level
    Ny = params.num_odes_y
    Nz = params.num_odes_z

    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    dy = L / Ny

    D_co2 = params.diffusion_coeff_co2_water    #m^2 / hr, depends on Temperature
    P_a = params.saturated_vapor_pressure_water_air(Tamb) 

    #density to M conversions
    molecular_weight_co2 = params.molecular_weight_co2*1000 #g/mol
 
    #initial concentrations


    Nelements = (Ny+1) * (Nz+1)
    Nelements1 = Ny+1
    C = max.(C, 1E-09) #don't allow negative values
    DIC = max.(DIC, 1E-09)
    dCO2 = zeros(Nelements, 1)
    dDIC = zeros(Nelements, 1)
    K(T,S) = params.dVavgdx(T,S,RH,WNDSPD,P_a)
    M_Evap(T) = params.evaporation_mass_flux(T,WNDSPD,RH,P_a)
    rho_solution(S) = params.density_solution(Tamb,S)
    S(T,x) = params.salinity(T,RH,WNDSPD,P_a,x)

   
    Hght(T,x,S) = params.height(T,S,RH,WNDSPD,P_a,x)
    dz(T,x,S) = Hght(T,x,S)/Nz

    ##Velocity Profile
    Re(H,T,V) = params.reynolds_number(H,T,V)
    V_profile = zeros(Nelements, 1)
    Vavg = zeros(Nelements1)

   
    dz_v = zeros(Ny + 1,1)
    Sal = zeros(Ny + 1, 1)
    Ht = zeros(Ny + 1,1)
    R_co2 = zeros(Nelements,1)
    mu_v = zeros(Nelements, 1)
    M = zeros(Nelements,1)


    M[pos2idx(0,0)] = W*dy*H*rho_solution(35)
    Vavg[pos2idx(0,0)] = params.volumetric_flow_rate_o/(W*H)

    for i in 1:Ny
        M[pos2idx(i,0)] = M[pos2idx(i-1,0)] - dy*(1/Vavg[pos2idx(i-1,0)])*M_Evap(Temperature[pos2idx(i-1,0)])*dy*W
        Vavg[pos2idx(i,0)] = sqrt((M[pos2idx(i-1,0)]*(Vavg[pos2idx(i-1,0)])^2)/(M[pos2idx(i,0)]))
    end

    for i in 0:Ny
 
                #salinity
            Sal[pos2idx(i,0)] = S(Temperature[pos2idx(i,0)],i) #kg/m3
                #height 
            #Vavg[pos2idx(i,0)] = params.avg_velocity(Temperature[pos2idx(i,0)], Sal[pos2idx(i,0)], RH, WNDSPD, P_a,i)

            Ht[pos2idx(i,0)] = Hght(Temperature[pos2idx(i,0)],i,Vavg[pos2idx(i,0)]) #m
                #increments in z direction
            dz_v[pos2idx(i,0)] = dz(Temperature[pos2idx(i,0)],i,Vavg[pos2idx(i,0)]) #m
    end

    co2_availability_factor(C_co2) = params.co2_availability_factor(C_co2)
   
   
    
    
    
    # ==========  DISSOLVED CARBON DIOXIDE BALANCE ========
    # C is dissolved CO2 [kg CO2 / m^3 water ]

    #BC1: CO2(y,z)  at y = 0 is constant [CO2_in]
    #     dCO2(y,z) at y = 0 is zero.

    wavelength = zeros(301,1)
    absorbance = zeros(301,1)
    percent_light = zeros(301,1)
    
    csv_reader1 = CSV.File("WAVE_ABS.csv", header=["col1", "col2","col3"])
    for (i,row) in enumerate(csv_reader1)
        wavelength[i] = convert(Float64, row.col1) #nm 400-700
        absorbance[i] = convert(Float64, row.col2) #m2/mol
        percent_light[i] = convert(Float64,row.col3)
    end

  

    watt_to_umolm2s = 0.425*4.6

    I_avg_vec = zeros(301,Nelements)

    for i = 1:301
        for j = 0:Ny
            I_avg_vec[i,pos2idx(j,0)] = GHI*watt_to_umolm2s*percent_light[i]
            for k = 1:Nz
                I_avg_vec[i,pos2idx(j,k)] =  I_avg_vec[i,pos2idx(j,k-1)]*exp(-absorbance[i]*(C_biomass[pos2idx(j,k)])*dz_v[pos2idx(j,0)])
                
            end
        end
    end

    I_avg = zeros(Nelements,1)
    Y_xphm = 1.38/24 #maximal biomass yield on light, molx/mol
    a_x = 4.81 #absorbance cross section, m2/mol

    phiL = zeros(Nelements,1)


    mw_co2 = params.molecular_weight_co2
    co2_to_M = (1/(mw_co2*1000)) #g/m3 to mol/m3
    pH = zeros(Nelements,1)

    phiCO2 = zeros(Nelements, 1)

    for i = 0:Ny
        for j = 0:Nz
            if  C[pos2idx(i,j)] <= 0.1782/0.0315 && C[pos2idx(i,j)] >= 0.689655
                phiCO2[pos2idx(i,j)] = (0.0261*C[pos2idx(i,j)] - 0.018)/0.129383
            elseif C[pos2idx(i,j)] > 0.1782/0.0315 && C[pos2idx(i,j)] <= 29.67
                phiCO2[pos2idx(i,j)] = (-0.0054*C[pos2idx(i,j)] + 0.1602)/0.129383
            else
                phiCO2[pos2idx(i,j)] = 0
            end
        end
    end

    for i = 0:Ny
        for j = 0:Nz
            I_avg[pos2idx(i,j)] = sum(I_avg_vec[1:301,pos2idx(i,j)])
            phiL[pos2idx(i,j)] = tanh((Y_xphm*a_x*I_avg[pos2idx(i,j)]*1E-06)/(bm*(1/3600))) #umol to mol
            pH[pos2idx(i,j)] = params.pH_interp(min(DIC[pos2idx(i,j)],70000))
            
            mu_v[pos2idx(i,j)] = phiCO2[pos2idx(i,j)]*phiL[pos2idx(i,j)]*mu_model(Temperature[pos2idx(i,j)],Sal[pos2idx(i,0)],C[pos2idx(i,j)])-0.003621  #1/hr
            #0.00361 is maintenance rate from Krishnan et al.
            
            R_co2[pos2idx(i,j)] = co2_per_biomass * mu_v[pos2idx(i,j)] * ((C_biomass[pos2idx(i,j)])*1000) #uptake of co2 by biomass, g/m3/hr

            if Re(Ht[pos2idx(0,0)],Temperature[pos2idx(0,0)],Vavg[pos2idx(0,0)]) > 2000
                V_profile[pos2idx(i,j)] = Vavg[pos2idx(i,0)]
            else
                V_profile[pos2idx(i,j)] = params.velocity_profile_lam(Vavg[pos2idx(i,0)],j,Ht[pos2idx(i,0)])
            end
            
            
        end
    end

    for i = 0:Ny
        C[pos2idx(i,0)] =  0.4401 #initial CO2
        dCO2[pos2idx(i,0)] = 0
        
    end

    #y = dy .. Ny,  z = 0 or z = Nz
    dMt(T,C_co2,H) = params.dMt(T,C_co2,H) 
    mf_co2 = params.mol_frac_co2
    dC_sparge = zeros(Nelements,1)
    G = params.G
    y_out(T,C,H) = params.y_out(T,C,H)
   

    mol_frac_co2 = zeros(Nelements,1)

    CO2_Conc = zeros(Nelements,1)

    for i = 0:Ny
        for j = 0:Nz
            CO2_Conc[pos2idx(i,j)] = C[pos2idx(i,j)]
        end
    end

    for i = 1:Ny
        mol_frac_co2[pos2idx(i,Nz)] = params.y_out(Temperature[pos2idx(i,Nz)],C[pos2idx(i,Nz)],dz_v[pos2idx(i,0)],mf_co2)
        dC_sparge[pos2idx(i,Nz)] = (G*(mf_co2 - mol_frac_co2[pos2idx(i,Nz)])*molecular_weight_co2)/(W*dy*dz_v[pos2idx(i,0)])
    end


    for i = 1:Ny
        for j = 1:Nz-1
            mol_frac_co2[pos2idx(i,Nz - j)] = params.y_out(Temperature[pos2idx(i,Nz - j)],C[pos2idx(i,Nz - j)],dz_v[pos2idx(i,0)],mol_frac_co2[pos2idx(i, (Nz - j) +1)])
            dC_sparge[pos2idx(i,Nz - j)] = (G*(mol_frac_co2[pos2idx(i, (Nz - j) +1)] - mol_frac_co2[pos2idx(i,Nz-j)])*molecular_weight_co2)/(W*dy*dz_v[pos2idx(i,0)])
        end
    end

    
    strip_pos = params.strip_position
    strip_Q = params.volumetric_flow_rate_strip #m3/hr
    strip_C = params.strip_concentration #g/m3

    Strip_add = zeros(Nelements,1)


    #all except for floor of reactor

    Q = zeros(Nelements,1)

    for i = 0:Ny
        for j = 0:Nz
            Q[pos2idx(i,j)] = V_profile[pos2idx(i,j)]*Ht[pos2idx(i,0)] #m2/hr
        end
    end
    
    
    for i=1:Ny
        for j=1:Nz

            Strip_add[pos2idx(strip_pos[1],j)] = (strip_Q*strip_C)/(Ht[pos2idx(strip_pos[1],0)]*dy*W)
            Strip_add[pos2idx(strip_pos[2],j)] = (strip_Q*strip_C)/(Ht[pos2idx(strip_pos[2],0)]*dy*W)
            Strip_add[pos2idx(strip_pos[3],j)] = (strip_Q*strip_C)/(Ht[pos2idx(strip_pos[3],0)]*dy*W)
            Strip_add[pos2idx(strip_pos[4],j)] = (strip_Q*strip_C)/(Ht[pos2idx(strip_pos[4],0)]*dy*W)
            Strip_add[pos2idx(strip_pos[5],j)] = (strip_Q*strip_C)/(Ht[pos2idx(strip_pos[5],0)]*dy*W)
            Strip_add[pos2idx(strip_pos[6],j)] = (strip_Q*strip_C)/(Ht[pos2idx(strip_pos[6],0)]*dy*W)

            dCO2[pos2idx(i,j)]=  (D_co2(Temperature[pos2idx(i,j)]) * (CO2_Conc[pos2idx(i-1,j)] + CO2_Conc[pos2idx(min(i+1,Ny),j)] - 2*CO2_Conc[pos2idx(i,j)]) / dy^2           #diffusion CO2 in y-direction
                                   + D_co2(Temperature[pos2idx(i,j)]) * (CO2_Conc[pos2idx(i,j-1)] + CO2_Conc[pos2idx(i,min(j+1,Nz))] - 2*CO2_Conc[pos2idx(i,j)]) / dz_v[pos2idx(i,0)]^2             #diffusion CO2 in z-direction
                                   - (Q[pos2idx(i,j)]*CO2_Conc[pos2idx(i,j)] - Q[pos2idx(i-1,j)]*CO2_Conc[pos2idx(i-1,j)]) / (dy*Ht[pos2idx(i,0)])                            #convection of CO2 in y-direction
                                   - (R_co2[pos2idx(i,j)])
                                   + Strip_add[pos2idx(i,j)])                                                 #consumption of CO2 by biomass growth
                                   #+ dC_sparge[pos2idx(i,j)])
                                                
        end
    end
    
    dDIC .= dCO2.*co2_to_M.*1000

    C_new = zeros(Nelements,1)

    
    for i = 0:Ny
        for j = 0:Nz
            DIC_new = max(DIC[pos2idx(i,j)]+ dDIC[pos2idx(i,j)],1E-09)
            C_new[pos2idx(i,j)] = params.CO2_interp(min(DIC_new,70000)) #calculate new CO2 after DIC dissociates
        end
    end
    
    
    for i = 1:Ny
        for j = 1:Nz
           dCO2[pos2idx(i,j)] = C_new[pos2idx(i,j)] - C[pos2idx(i,j)] #this sets the new CO2 value to C_new
        end
    end
    
   
    @views dX[1+2*Nelements:3*Nelements] .= dCO2        # order matters! The @views operator takes a slice out of an array without making a copy.
    @views dX[1+3*Nelements:4*Nelements] .= dDIC
    nothing
end




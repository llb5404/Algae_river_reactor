

function PDE_AlgaeBiomass!(dX, C, DIC,CO2,Temperature, params, t)
    #Units for params can be found in LoadParameters.jl
    
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data

    t_hour1 = floor(Int64, t)
    t_hour2 = floor(Int64, t)+1

    CO2 = max.(CO2, 1E-09)
    DIC = max.(DIC, 1E-09)
    
    data_begin = params.data_begin

    GHI = GHI_Data[data_begin + t_hour1] * (t-t_hour1) + GHI_Data[data_begin + t_hour2] * (t_hour2 - t)
    Tamb = max.(((Tamb_Data[data_begin + t_hour1] * (t-t_hour1) + Tamb_Data[data_begin + t_hour2] * (t_hour2 - t))),0)
    RH = RH_Data[data_begin + t_hour1] * (t-t_hour1) + RH_Data[data_begin + t_hour2] * (t_hour2 - t)
    WNDSPD = max.((WNDSPD_Data[data_begin + t_hour1] * (t-t_hour1) + WNDSPD_Data[data_begin + t_hour2] * (t_hour2 - t)),0)
    
    Dy = params.biomass_diffusion_coefficient_y
    Dz = params.biomass_diffusion_coefficient_z
    L = params.reactor_length 
    H = params.reactor_initial_liquid_level 
    W = params.reactor_width
    Ny = params.num_odes_y 
    Nz = params.num_odes_z
    dy = L / Ny
 


    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]


    Nelements = (Ny+1) * (Nz+1)
    Nelements1 = Ny + 1
    
    C = max.(C, 1E-09) #don't allow negative values
    dC = zeros( Nelements, 1)
    P_a = params.saturated_vapor_pressure_water_air(Tamb) 
    K(T,S) = params.dVavgdx(T,S,RH,WNDSPD,P_a)
    M_Evap(T) = params.evaporation_mass_flux(T,WNDSPD,RH,P_a)
    rho_solution(S) = params.density_solution(Tamb,S)
    S(T,x) = params.salinity(T,RH,WNDSPD,P_a,x)

    co2_availability_factor(Co2) = params.co2_availability_factor(Co2)
    temperature_factor(T) = params.temperature_factor(T)
    salinity_factor(S) = params.salinity_factor(S)
    bm = params.max_biomass_specific_growth_rate
    
    mu(T,S,Co2) = params.biomass_specific_growth_rate(T, S, Co2)
    

    Hght(T,x,S) = params.height(T,S,RH,WNDSPD,P_a,x)
    dz(T,x,S) = Hght(T,x,S)/Nz

   # Velocity Profile
    Re(H,T,V) = params.reynolds_number(H,T,V)
    V_profile = zeros(Nelements, 1)
    dz_v = zeros(Nelements,1)
    Sal = zeros(Nelements, 1)
    Ht = zeros(Nelements,1)
    Vavg = zeros(Nelements1)

    M = zeros(Nelements,1)


    M[pos2idx(0,0)] = W*dy*H*rho_solution(35)
    Vavg[pos2idx(0,0)] = params.volumetric_flow_rate_o/(W*H)

    for i in 1:Ny
        M[pos2idx(i,0)] = M[pos2idx(i-1,0)] - dy*(1/Vavg[pos2idx(i-1,0)])*M_Evap(Temperature[pos2idx(i-1,0)])*dy*W
        Vavg[pos2idx(i,0)] = (M[pos2idx(i-1,0)]*Vavg[pos2idx(i-1,0)])/(M[pos2idx(i,0)])
    end


    for i in 0:Ny
        for j in 0:Nz

                #salinity
                Sal[pos2idx(i,j)] = S(Temperature[pos2idx(i,j)],i) #kg/m3
                #height 

                #Vavg[pos2idx(i,0)] = params.avg_velocity(Temperature[pos2idx(i,0)], Sal[pos2idx(i,0)], RH, WNDSPD, P_a,i)

                Ht[pos2idx(i,j)] = Hght(Temperature[pos2idx(i,j)],i,Vavg[pos2idx(i,0)]) #m
                #increments in z direction
                dz_v[pos2idx(i,j)] = dz(Temperature[pos2idx(i,j)],i,Vavg[pos2idx(i,0)]) #m

                if Re(Ht[pos2idx(0,0)],Temperature[pos2idx(0,0)],Vavg[pos2idx(0,0)]) > 2000
                    V_profile[pos2idx(i,j)] = Vavg[pos2idx(i,0)]
                else
                    V_profile[pos2idx(i,j)] = params.velocity_profile_lam(Vavg[pos2idx(i,0)],j,Ht[pos2idx(i,j)])
                end
         
        end
    end

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
                I_avg_vec[i,pos2idx(j,k)] = I_avg_vec[i,pos2idx(j,k-1)]*exp(-absorbance[i]*(C[pos2idx(j,k)])*dz_v[pos2idx(j,k)]) #the density of the solution (not mass, has an effect on light)
            end
        end
    end

    I_avg = zeros(Nelements,1)

    for i = 0:Ny
        for j = 0:Nz
            I_avg[pos2idx(i,j)] = sum(I_avg_vec[1:301,pos2idx(i,j)])
        
        end
    end

    Y_xphm = 1.38/24 #maximal biomass yield on light, molx/mol
    a_x = 4.81 #absorbance cross section, m2/mol

    phiL = zeros(Nelements,1)

    
    for i = 0:Ny
        for j = 0:Nz
            phiL[pos2idx(i,j)] = tanh((Y_xphm*a_x*I_avg[pos2idx(i,j)]*1E-06)/(bm*(1/3600))) #umol to mol
        end
    end

    phiCO2 = zeros(Nelements, 1)

    for i = 0:Ny
        for j = 0:Nz
            if  CO2[pos2idx(i,j)] <= 0.1782/0.0315 && CO2[pos2idx(i,j)] >= 0.689655
                phiCO2[pos2idx(i,j)] = (0.0261*CO2[pos2idx(i,j)] - 0.018)/0.129383
            elseif CO2[pos2idx(i,j)] > 0.1782/0.0315 && CO2[pos2idx(i,j)] <= 29.67
                phiCO2[pos2idx(i,j)] = (-0.0054*CO2[pos2idx(i,j)] + 0.1602)/0.129383
            else
                phiCO2[pos2idx(i,j)] = 0
            end
        end
    end

    

    
    mw_co2 = params.molecular_weight_co2
    co2_to_M = (1/(mw_co2*1000)) #g/m3 to mol/m3

    #Vector of specific growth rate values
    mu_v = zeros(Nelements, 1)
    pH = zeros(Nelements,1)
    for i in 0:Ny
        for j in 0:Nz

            pH[pos2idx(i,j)] = params.pH_interp(min(DIC[pos2idx(i,j)],70000))

            mu_v[pos2idx(i,j)] = phiL[pos2idx(i,j)]*phiCO2[pos2idx(i,j)]*mu(Temperature[pos2idx(i,j)],Sal[pos2idx(i,j)],CO2[pos2idx(i,j)]) -0.003621
    
            
        end
    end
    



   
    


    # ==========  ALGAE GROWTH ========

    #BC1: C(y,z) at y = 0 is constant [Cin]
    #     dC(y,z) at y = 0 is zero.

    #y = 0, z = 0 ... Nz

    #BC1: C(y,z) at y = 0 is constant [Cin]
    #     dC(y,z) at y = 0 is zero.
    #y = 0, z = 0 ... Nz
    dC[pos2idx(0,0:Nz)] .= 0.0
    #BC2: molar flux is zero at z = 0, H
    # C[pos2idx(i,-1)] = C[pos2idx(i,0)]
    # C[pos2idx(i,Nz+1)] = C[pos2idx(i,Nz)]
    # C[pos2idx(0,-1)] = C[pos2idx(i,0)]
    # C[pos2idx(0,Nz+1)] = C[pos2idx(i,Nz)]
    # C[pos2idx(Ny+1,z)] = C[pos2idx(Ny,z)]
    #BC3: at y = Ny, convection only [accumulation of biomass at end]
    #y = Ny,  z = 0 .. Nz
    vT(T,S) = params.vT(T,S)
    Q = zeros(Nelements,1)

    for i = 0:Ny
        for j = 0:Nz
            Q[pos2idx(i,j)] = V_profile[pos2idx(i,j)]*Ht[pos2idx(i,0)] #m2/hr
        end
    end

    C_conc = zeros(Nelements,1)



    for i = 0:Ny
        for j = 0:Nz
            C_conc[pos2idx(i,j)] = C[pos2idx(i,j)]
        end
    end

    #dC/dy: C[pos2idx(i+1,j)] = 
    
    for j=0:Nz
        

        dC[pos2idx(Ny,j)] = ( (+ Dz * (C_conc[pos2idx(Ny,max(0,j-1))] + C_conc[pos2idx(Ny,min(Nz,j+1))] - 2*C_conc[pos2idx(Ny,j)]) / dz_v[pos2idx(Ny,j)]^2)
                            - (Q[pos2idx(Ny,j)]*C_conc[pos2idx(Ny,j)] - Q[pos2idx(Ny-1,j)]*C_conc[pos2idx(Ny-1,j)]) / (dy*Ht[pos2idx(Ny,j)]) #change in concentration due to height change
                            )
        
    end

    #y = dy .. Ny-dy,  z = 0 or z = Nz
    for i=1:Ny-1
        

        dC[pos2idx(i,0)] =     ( (Dy * (C_conc[pos2idx(i-1,0)] + C_conc[pos2idx(i+1,0)] - 2*C_conc[pos2idx(i,0)]) / dy^2  #small
                               + Dz * (C_conc[pos2idx(i,0)] + C_conc[pos2idx(i,0+1)] - 2*C_conc[pos2idx(i,0)]) /dz_v[pos2idx(i,0)]^2)
                               - (Q[pos2idx(i,0)]*C_conc[pos2idx(i,0)] - Q[pos2idx(i-1,0)]*C_conc[pos2idx(i-1,0)]) / (dy*Ht[pos2idx(i,0)])
                               + mu_v[pos2idx(i,0)] * C_conc[pos2idx(i,0)]
                               )


        dC[pos2idx(i,Nz)] =    ( (Dy * (C_conc[pos2idx(i-1,Nz)] + C_conc[pos2idx(i+1,Nz)] - 2*C_conc[pos2idx(i,Nz)]) / dy^2 #small
                               + Dz * (C_conc[pos2idx(i,Nz-1)] + C_conc[pos2idx(i,Nz)] - 2*C_conc[pos2idx(i,Nz)]) / dz_v[pos2idx(i,Nz)]^2)
                               - (Q[pos2idx(i,Nz)]*C_conc[pos2idx(i,Nz)] - Q[pos2idx(i-1,Nz)]*C_conc[pos2idx(i-1,Nz)]) / (dy*Ht[pos2idx(i,Nz)])
                               + mu_v[pos2idx(i,Nz)] * C_conc[pos2idx(i,Nz)]
                               )
    end

    for i=1:Ny-1
        for j=1:Nz-1

            dC[pos2idx(i,j)] = ( (Dy * (C_conc[pos2idx(i-1,j)] + C_conc[pos2idx(i+1,j)] - 2*C_conc[pos2idx(i,j)]) / dy^2
                               + Dz * (C_conc[pos2idx(i,j-1)] + C_conc[pos2idx(i,j+1)] - 2*C_conc[pos2idx(i,j)]) / dz_v[pos2idx(i,j)]^2)
                               - (Q[pos2idx(i,j)]*C_conc[pos2idx(i,j)] - Q[pos2idx(i-1,j)]*C_conc[pos2idx(i-1,j)]) / (dy*Ht[pos2idx(i,j)])
                               + mu_v[pos2idx(i,j)] * C_conc[pos2idx(i,j)]
                               )
                  
        end
    end

    C = max.(C, 1E-09)
    
    @views dX[1:Nelements] .= dC        # order matters! The @views operator takes a slice out of an array without making a copy.

    nothing
end
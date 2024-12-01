

function PDE_AlgaeBiomass!(dX, C, DIC,CO2,Temperature, params, t)
    #Units for params can be found in LoadParameters.jl
    
    #Environmental Data
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data

    #Timestep
    t_hour1 = floor(Int64, t)
    t_hour2 = floor(Int64, t)+1

    #Ensures that CO2 and DIC do not reach 0
    CO2 = max.(CO2, 1E-09)
    DIC = max.(DIC, 1E-09)
    
    data_begin = params.data_begin

    #Adjusts Environmental Data for each timestep
    GHI = GHI_Data[data_begin + t_hour1] * (t-t_hour1) + GHI_Data[data_begin + t_hour2] * (t_hour2 - t)
    Tamb = max.(((Tamb_Data[data_begin + t_hour1] * (t-t_hour1) + Tamb_Data[data_begin + t_hour2] * (t_hour2 - t))),0)
    RH = RH_Data[data_begin + t_hour1] * (t-t_hour1) + RH_Data[data_begin + t_hour2] * (t_hour2 - t)
    WNDSPD = max.((WNDSPD_Data[data_begin + t_hour1] * (t-t_hour1) + WNDSPD_Data[data_begin + t_hour2] * (t_hour2 - t)),0)
    
    #Biomass Diffusion
    Dy = params.biomass_diffusion_coefficient_y
    Dz = params.biomass_diffusion_coefficient_z

    #Geometric Properties
    L = params.reactor_length 
    H = params.reactor_initial_liquid_level 
    W = params.reactor_width
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
    

    C = max.(C, 1E-09) #biomass (g/L) don't allow negative values
    dC = zeros(Nelements, 1) #dC/dt

    #Other parameters
    P_a = params.saturated_vapor_pressure_water_air(Tamb) 
    M_Evap = params.evaporation_mass_flux(Tamb,WNDSPD,RH,P_a)
    rho_water = params.density_water(Tamb)
    S(x) = params.salinity(Tamb,RH,WNDSPD,P_a,x)
    bm = params.max_biomass_specific_growth_rate
    mu(T,S) = params.biomass_specific_growth_rate(T, S)
    Re(H,T,V) = params.reynolds_number(H,T,V)

    #Initialize vectors for vertical increments, salinity, height
    dz_v = zeros(Nelements1)
    Sal = zeros(Nelements1)
    Ht = zeros(Nelements1)

    #Vavg determined assuming conservation of kinetic energy (evaporated water has 0 kinetic energy)
    M = zeros(Nelements1) #mass of water in vertical slice (kg)
    Vavg = zeros(Nelements1) #average velocity

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
                if Re(Ht[pos2idx(0,0)],Temperature[pos2idx(0,0)],Vavg[pos2idx(0,0)]) > 2000 #turbulent
                    V_profile[pos2idx(i,j)] = Vavg[pos2idx(i,0)]
                else #laminar
                    V_profile[pos2idx(i,j)] = params.velocity_profile_lam(Vavg[pos2idx(i,0)],j,Ht[pos2idx(i,0)])
                end
         
        end
    end

    #volumetric flowrate and concentration adjustment factor (used to adjust biomass and co2 concentrations due to added water, affects light penetration and CO2 concentration terms in mu)
    Q = zeros(Nelements,1) #volumetric flowrate (m3/hr)
    Q_adj = zeros(Nelements,1) #concentration adjustment factor

    for i = 0:Ny
        for j = 0:Nz
            Q[pos2idx(i,j)] = Vavg[pos2idx(i,0)]*Ht[pos2idx(i,0)] #m2/hr
            Q_adj[pos2idx(i,j)] = Q[pos2idx(0,j)]/Q[pos2idx(i,j)]
        end
    end

    #Calculation of light adjustment factor for biomass growth
    wavelength = zeros(301,1)
    absorbance = zeros(301,1)
    percent_light = zeros(301,1)
    
    csv_reader1 = CSV.File("WAVE_ABS.csv", header=["col1", "col2","col3"])
    for (i,row) in enumerate(csv_reader1)
        wavelength[i] = convert(Float64, row.col1) #wavelength of light, nm 400-700
        absorbance[i] = convert(Float64, row.col2) #absorbance of wavelength of light by P celeri, m2/mol
        percent_light[i] = convert(Float64,row.col3) #percentage of sunlight in each wavelength
    end

    watt_to_umolm2s = 0.425*4.6 #watt to umolm2s conversion

    #average light intensity for each wavelength
    I_avg_vec = zeros(301,Nelements)

    for i = 1:301
        for j = 0:Ny
            I_avg_vec[i,pos2idx(j,0)] = GHI*watt_to_umolm2s*percent_light[i]
            for k = 1:Nz
                I_avg_vec[i,pos2idx(j,k)] = I_avg_vec[i,pos2idx(j,k-1)]*exp(-absorbance[i]*(Q_adj[pos2idx(j,k)]*C[pos2idx(j,k)])*dz_v[pos2idx(j,0)]) 
            end
        end
    end

    #average light intensity summed across all wavelengths
    I_avg = zeros(Nelements,1)

    for i = 0:Ny
        for j = 0:Nz
            I_avg[pos2idx(i,j)] = sum(I_avg_vec[1:301,pos2idx(i,j)])
        
        end
    end

    Y_xphm = 1.38/24 #maximal biomass yield on light, molx/mol
    a_x = 4.81 #absorbance cross section, m2/mol

    #Growth adjustment factor - light
    phiL = zeros(Nelements,1)

    
    for i = 0:Ny
        for j = 0:Nz
            phiL[pos2idx(i,j)] = tanh((Y_xphm*a_x*I_avg[pos2idx(i,j)]*1E-06)/(bm*(1/3600)))
        end
    end


    #Growth adjustment factor - CO2
    phiCO2 = zeros(Nelements, 1)

    #Piecewise linear function obtained from experiments on P celeri (see documentation)
    for i = 0:Ny
        for j = 0:Nz
            if  CO2[pos2idx(i,j)]*Q_adj[pos2idx(i,j)] <= 0.1782/0.0315 && CO2[pos2idx(i,j)]*Q_adj[pos2idx(i,j)] >= 0.689655
                phiCO2[pos2idx(i,j)] = (0.0261*CO2[pos2idx(i,j)]*Q_adj[pos2idx(i,j)] - 0.018)/0.129383
            elseif CO2[pos2idx(i,j)]*Q_adj[pos2idx(i,j)] > 0.1782/0.0315 && CO2[pos2idx(i,j)]*Q_adj[pos2idx(i,j)] <= 29.67
                phiCO2[pos2idx(i,j)] = (-0.0054*CO2[pos2idx(i,j)]*Q_adj[pos2idx(i,j)] + 0.1602)/0.129383
            else
                phiCO2[pos2idx(i,j)] = 0
            end
        end
    end

    #Biomass specific growth term
    mu_v = zeros(Nelements, 1)
    for i in 0:Ny
        for j in 0:Nz

            mu_v[pos2idx(i,j)] = phiCO2[pos2idx(i,j)]*phiL[pos2idx(i,j)]*mu(Temperature[pos2idx(i,j)],Sal[pos2idx(i,0)]) -0.003621 

        end
    end
    



   
    


    # ==========  ALGAE GROWTH ========

    #BC1: C(y,z) at y = 0 is constant [Cin]
    dC[pos2idx(0,0:Nz)] .= 0.0

    #BC 2: No biomass growth (mu) at outlet of reactor (Ny)
    for j=0:Nz
        

        dC[pos2idx(Ny,j)] = ( (+ Dz * (C[pos2idx(Ny,max(0,j-1))] + C[pos2idx(Ny,min(Nz,j+1))] - 2*C[pos2idx(Ny,j)]) / dz_v[pos2idx(Ny,0)]^2)
                            - (V_profile[pos2idx(Ny,j)]*(C[pos2idx(Ny,j)] - C[pos2idx(Ny-1,j)])/dy) #change in concentration due vol flr change
                            )
        
    end

    #BC 3: At floor and surface, no diffusion of biomass in the vertical direction 
    for i=1:Ny-1
        

        dC[pos2idx(i,0)] = ( (Dy * (C[pos2idx(i-1,0)] + C[pos2idx(i+1,0)] - 2*C[pos2idx(i,0)]) / dy^2  #small
        + Dz * (C[pos2idx(i,0)] + C[pos2idx(i,0+1)] - 2*C[pos2idx(i,0)]) /dz_v[pos2idx(i,0)]^2)
        - (V_profile[pos2idx(i,0)]*(C[pos2idx(i,0)] - C[pos2idx(i-1,0)])/dy) 
        + mu_v[pos2idx(i,0)] * C[pos2idx(i,0)]
        )


        dC[pos2idx(i,Nz)] =    ( (Dy * (C[pos2idx(i-1,Nz)] + C[pos2idx(i+1,Nz)] - 2*C[pos2idx(i,Nz)]) / dy^2 #small #diffusion in y direction changes slightly
                               + Dz * (C[pos2idx(i,Nz-1)] + C[pos2idx(i,Nz)] - 2*C[pos2idx(i,Nz)]) / dz_v[pos2idx(i,0)]^2)
                               - (V_profile[pos2idx(i,Nz)]*(C[pos2idx(i,Nz)] - C[pos2idx(i-1,Nz)])/dy) #rate of mass transfer due to velocity is also the same
                               + mu_v[pos2idx(i,Nz)] * C[pos2idx(i,Nz)] #this remains because rate of mass increase is the same regardless of volume
                               )
    end

    for i=1:Ny-1
        for j=1:Nz-1

            dC[pos2idx(i,j)] = ( (Dy * (C[pos2idx(i-1,j)] + C[pos2idx(i+1,j)] - 2*C[pos2idx(i,j)]) / dy^2 
                               + Dz * (C[pos2idx(i,j-1)] + C[pos2idx(i,j+1)] - 2*C[pos2idx(i,j)]) / dz_v[pos2idx(i,0)]^2)
                               - (V_profile[pos2idx(i,j)]*(C[pos2idx(i,j)] - C[pos2idx(i-1,j)])/dy)
                               + mu_v[pos2idx(i,j)] * C[pos2idx(i,j)]
                               )
                  
        end
    end

    C = max.(C, 1E-09)
    
    @views dX[1:Nelements] .= dC        # order matters! The @views operator takes a slice out of an array without making a copy.

    nothing
end
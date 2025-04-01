
function PDE_CO2!(dX, C, DIC, C_biomass, Temperature, A, params, t)
    #Units for params can be found in LoadParameters.jl

    #Environmental Data
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data

    #Timestep
    t_hour1 = floor(Int64, t)
    t_hour2 = floor(Int64, t)+1

    data_begin = params.data_begin

    #Adjusts Environmental Data for each timestep
    GHI = GHI_Data[data_begin + t_hour1] * (t-t_hour1) + GHI_Data[data_begin + t_hour2] * (t_hour2 - t)
    Tamb = max.(((Tamb_Data[data_begin + t_hour1] * (t-t_hour1) + Tamb_Data[data_begin + t_hour2] * (t_hour2 - t))),0)
    RH = RH_Data[data_begin + t_hour1] * (t-t_hour1) + RH_Data[data_begin + t_hour2] * (t_hour2 - t)
    WNDSPD = max.((WNDSPD_Data[data_begin + t_hour1] * (t-t_hour1) + WNDSPD_Data[data_begin + t_hour2] * (t_hour2 - t)),0)

    mu(T,S) = params.biomass_specific_growth_rate(T, S)
    co2_per_biomass = params.co2_per_biomass
    bm = params.max_biomass_specific_growth_rate

    #Geometric Properties
    L = params.reactor_length
    W = params.reactor_width
    Nx = params.num_odes_x
    Ny = params.num_odes_y
    dx = L / Nx

    #Reindexing of 2D position in reactor to be contained in 1D vector
    pos2idx(y) = (y.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Nx+1) * floor( (pos-1) ./ (Nx.+1))), Integer(floor( (pos-1) ./ (Nx.+1)))]

    #Number of total elements in vector
    Nelements = (Nx+1)
  
    C = max.(C, 1E-09) #CO2 concentration (g/m3), don't allow 0 and negative values
    DIC = max.(DIC, 1E-09) #DIC concentration (g/m3), don't allow 0 and negative values
    dCO2 = zeros(Nelements, 1) #dCO2/dt
    dDIC = zeros(Nelements, 1) #dDIC/dt
    dA = zeros(Nelements)


    #Other parameters
    D_co2 = params.diffusion_coeff_co2_water
    P_a = params.saturated_vapor_pressure_water_air(Tamb) 
    S(Area) = params.salinity(Area)
    alpha = params.alpha
    m = params.m

    #Initialize vectors for vertical increments, salinity, height

    Sal = zeros(Nelements, 1)
    Vavg = zeros(Nelements)
    H = zeros(Nelements)
   
    for i in 0:Nx
        Vavg[pos2idx(i)] = alpha*A[pos2idx(i)]^(m-1)
        H[pos2idx(i)] = A[pos2idx(i)]/(W)
        #salinity
        Sal[pos2idx(i)] = S(A[pos2idx(i)]) #kg/m3
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
    I_avg = zeros(Nelements)

    for i = 1:301
        for j = 0:Nx
            I_avg_vec[i,pos2idx(j)] = GHI*watt_to_umolm2s*percent_light[i]*exp(-absorbance[i]*(C_biomass[pos2idx(j)])*H[pos2idx(j)]/2) 
        end
    end

    for i = 0:Nx
            I_avg[pos2idx(i)] = sum(I_avg_vec[1:301,pos2idx(i)])
    end

    Y_xphm = 1.38/24 #maximal biomass yield on light, molx/mol
    a_x = 4.81 #absorbance cross section, m2/mol

    #Growth adjustment factor - light
    phiL = zeros(Nelements,1)

    
    for i = 0:Nx
        phiL[pos2idx(i)] = tanh((Y_xphm*a_x*I_avg[pos2idx(i)]*1E-06)/(bm*(1/3600)))
    end


    #Growth adjustment factor - CO2
    phiCO2 = zeros(Nelements, 1)

    #Piecewise linear function obtained from experiments on P celeri (see documentation)
    for i = 0:Nx
            if  C[pos2idx(i)] <= 0.1782/0.0315 && C[pos2idx(i)] >= 0.689655
                phiCO2[pos2idx(i)] = (0.0261*C[pos2idx(i)] - 0.018)/0.129383
            elseif C[pos2idx(i)] > 0.1782/0.0315 && C[pos2idx(i)] <= 29.67
                phiCO2[pos2idx(i)] = (-0.0054*C[pos2idx(i)] + 0.1602)/0.129383
            else
                phiCO2[pos2idx(i)] = 0
            end
    end

    #Biomass specific growth term (1/hr)
    mu_v = zeros(Nelements, 1)
    #CO2 uptake term (g/m3/hr)
    R_co2 = zeros(Nelements,1)
    for i in 0:Nx
            mu_v[pos2idx(i)] = phiCO2[pos2idx(i)]*phiL[pos2idx(i)]*mu(Temperature[pos2idx(i)],Sal[pos2idx(i)]) -0.003621 
            R_co2[pos2idx(i)] = max(co2_per_biomass * mu_v[pos2idx(i)] * ((C_biomass[pos2idx(i)])*1000),0) #g/m3/hr (when algae die, they release CO2)
    end

    
    CO2_star(T) = params.C_CO2_star(T)
    k_co2(T) = params.k_co2(WNDSPD,T)

    # ==========  DISSOLVED CARBON DIOXIDE BALANCE ========
    # C is dissolved CO2 [kg CO2 / m^3 water ]

    #BC1: CO2(y,z)  at y = 0 is constant [CO2_in]
    #     dCO2(y,z) at y = 0 is zero.

    C_o = params.co2_init

    C[pos2idx(0)] =  C_o
   
    strip_Q(phiL) = params.volumetric_flow_rate_strip(phiL) #m3/hr
    strip_C = params.strip_concentration #g/m3
    q(Temp,phiL) = params.lat_flow(Temp,WNDSPD,RH,P_a,phiL)

    S_conc_new = zeros(Nelements)
    CO2_uptake = zeros(Nelements)
    CO2_str = zeros(Nelements)
    CO2_conv = zeros(Nelements)
    Term_adj = zeros(Nelements)
    CO2_flux = zeros(Nelements)

    for i=0:Nx
        dA[pos2idx(i)] = q(Temperature[pos2idx(i)],I_avg[pos2idx(i)]) -alpha*m*(((A[pos2idx(i)] + A[pos2idx(max(i-1,0))])/2)^(m-1))*(A[pos2idx(i)]-A[pos2idx(max(i-1,0))])/dx
    end
 
    beta = 100

    for i=1:Nx

            CO2_uptake[pos2idx(i)] = -R_co2[pos2idx(i)]
            CO2_conv[pos2idx(i)] = -(A[pos2idx(i)]*Vavg[pos2idx(i)]*C[pos2idx(i)] - A[pos2idx(i-1)]*Vavg[pos2idx(i-1)]*C[pos2idx(i-1)])/(A[pos2idx(i)]*dx)
            Term_adj[pos2idx(i)] = -dA[pos2idx(i)]*C[pos2idx(i)]/A[pos2idx(i)]
            CO2_flux[pos2idx(i)] = k_co2(Temperature[pos2idx(i)])*(CO2_star(Temperature[pos2idx(i)]) - C[pos2idx(i)])/H[pos2idx(i)]
            S_conc_new[pos2idx(i)] = 1500
            CO2_str[pos2idx(i)] = (strip_Q(I_avg[pos2idx(i)])*S_conc_new[pos2idx(i)])/(A[pos2idx(i)]*L)
        
            dCO2[pos2idx(i)]=  (+ CO2_uptake[pos2idx(i)]
                                   + CO2_str[pos2idx(i)]
                                   + CO2_conv[pos2idx(i)]
                                   + Term_adj[pos2idx(i)]
                                   + CO2_flux[pos2idx(i)])
    end
    
    #change in total dissolved carbon (umol) equal to change of CO2
    
    mw_co2 = 44.01 #g/mol
    co2_to_M = (1/(mw_co2*1000)) #g/m3 to mol/L
    M_to_uM = 1E6 # mol/L to umol/L
    

    dDIC .= dCO2.*co2_to_M.*M_to_uM #converts g/m3 to umol/L

    C_new = zeros(Nelements,1)
    for i = 0:Nx
       
            DIC_new = max(DIC[pos2idx(i)]+ dDIC[pos2idx(i)],1E-09) #new DIC concentration
            C_new[pos2idx(i)] = params.CO2_interp(min(DIC_new,70000)) #calculate new CO2 after DIC dissociates (seawater buffer)

    end
    
    for i = 1:Nx
      
        dCO2[pos2idx(i)] = C_new[pos2idx(i)] - C[pos2idx(i)] #this sets the new CO2 value to C_new
      
    end
   
    @views dX[1+2*Nelements:3*Nelements] .= dCO2        # order matters! The @views operator takes a slice out of an array without making a copy.
    @views dX[1+3*Nelements:4*Nelements] .= dDIC
    nothing
end




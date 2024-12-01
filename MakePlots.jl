
using CSV, Tables

function Plot_Biomass_Profile(Mout, CO2_out, DIC_out,Tout, T, params)


    Ny = params.num_odes_y
    Nz = params.num_odes_z
    L = params.reactor_length
    RL = params.real_length #Corresponds to 'actual length' of river reactor (outlet value taken at this position)
    Ny_new = trunc(Int,(RL/L)*Ny) #Index corresponding to 'actual length'
    dy = L/Ny
    W = params.reactor_width                    
    H = params.reactor_initial_liquid_level     
    TL = length(T) #number of elements in timestep

    #Number of total elements in vector
    Nelements = (Ny+1) * (Nz+1) # accounts for length and height
    NelementsT2 = (Nelements)*(TL+1) #accounts for length, height, time
    Nelements1 = Ny + 1 #accounts for length

    #Reindexing of X-Dim position in reactor to be contained in 1D vector
    pos2idx1(z) = z.+1 #1-D to 1-D
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1) #2-D to 1-D
    tpos2idx2(t,y,z) = Nelements*t + pos2idx(y,z) #3-D to 1-D
    
    #Environmental Data
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data
    data_begin = params.data_begin
 
    GHIout = zeros(TL,1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        GHIout[i] = GHI_Data[data_begin + t_hour1] * (T[i]-t_hour1) + GHI_Data[data_begin + t_hour2] * (t_hour2 - T[i])
    end
    WNDSPDout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        WNDSPDout[i] = WNDSPD_Data[data_begin + t_hour1]*(T[i]-t_hour1)+WNDSPD_Data[data_begin + t_hour2]*(t_hour2-T[i])
    end
    RHout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        RHout[i] = RH_Data[data_begin + t_hour1]*(T[i]-t_hour1)+RH_Data[data_begin + t_hour2]*(t_hour2-T[i])
    end
    Tambout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        Tambout[i] = (Tamb_Data[data_begin + t_hour1]*(T[i]-t_hour1)+Tamb_Data[data_begin + t_hour2]*(t_hour2-T[i]))
    end
    T_days = T ./ 24.0
    
    P_a(Tamb) = params.saturated_vapor_pressure_water_air(Tamb) #Pa
    Pa_out = zeros(TL, 1) #Pa
    for i in 1:TL
        Pa_out[i] = P_a(Tambout[i]) #pressure at ambient temp, Pa
    end

    ## Evaporation Loss Calculations

    Evap_C = zeros(TL, Ny+1) #evaporation constant, dimless
    X_air = zeros(TL, Ny+1) #ambient humidity ratio, dimless kg H2O/kg dry air
    P_atm = params.reference_pressure #Pa

    for i in 1:TL
        Evap_C[i,pos2idx(0:Ny,0)] .= 25 + 19*WNDSPDout[i]
        X_air[i,pos2idx(0:Ny,0)] .= 0.62198*(Pa_out[i]*(RHout[i]/100))/(P_atm-(Pa_out[i]*(RHout[i]/100)))
    end
   
    X_surface = zeros(TL, Ny+1) #surface humidity ratio, dimless
    M_Evap = zeros(TL, Ny+1) #evaporation mass flux, kg H2O/m2/hr
    for i in 1:TL
        for j in 0:Ny
           
            X_surface[i,pos2idx(j,0)] = 0.62198*(Pa_out[i]/(P_atm-Pa_out[i])) #kg H2O/kg dry air
            M_Evap[i,pos2idx(j,0)] = Evap_C[i,pos2idx(j,0)]*(X_surface[i,pos2idx(j,0)]-X_air[i, pos2idx(j,0)]) #kg H2O/m2-hr
        end
    end

    Evap_Loss_T = zeros(TL+1) #loss of water due to evaporation for each time increment
    for i = 1:TL
        Evap_Loss_T[i] = W*RL*M_Evap[i,pos2idx(Ny_new,0)]*(T[i]-T[max(i-1, 1)]) #selection of j index does not matter, M_evap is constant across reactor
    end

    Cum_Evap_Loss = sum(Evap_Loss_T) #sums across all time increments

    #Functions for height, vertical increment, volumetric flowrate, reynolds number, salinity
    Hght(T,S,RH,WNDSPD,P,x) = params.height(T,S,RH,WNDSPD,P,x)
    dz(T,S,RH,WNDSPD,P,x) = Hght(T,S,RH,WNDSPD,P,x)/Nz
    Q(T,WNDSPD,RH,P,x) = params.volumetric_flow_rate(T,WNDSPD,RH,P,x)
    Re(H,T,V) = params.reynolds_number(H,T,V)
    Sal(T,RH, WNDSPD,P,x) = params.salinity(T,RH,WNDSPD,P,x)

    #Initialize vectors
    Sout = zeros(TL, Ny+1)
    Ht = zeros(TL,Ny+1)
    dz_v = zeros(TL,Ny+1)
    V_prof_out = zeros(TL,Nelements)

    #Vavg determined assuming conservation of kinetic energy (evaporated water has 0 kinetic energy)
    Vavg = zeros(TL, Nelements1) #mass of water in vertical slice (kg)
    M = zeros(TL, Nelements1) #average velocity

    for i in 1:TL
        M[i, pos2idx(0,0)] = W*dy*H*params.density_water(Tout[i,pos2idx(0,0)])
        Vavg[i, pos2idx(0,0)] = params.volumetric_flow_rate_o/(W*H)
    end

    for i in 1:TL
        for j in 1:Ny
            M[i, pos2idx(j,0)] = M[i, pos2idx(j-1,0)] - dy*(1/Vavg[i, pos2idx(j-1,0)])*params.evaporation_mass_flux(Tout[i,pos2idx(j-1,0)],WNDSPDout[i],RHout[i],Pa_out[i])*dy*W
            Vavg[i, pos2idx(j,0)] = sqrt((M[i, pos2idx(j-1,0)]*(Vavg[i, pos2idx(j-1,0)])^2)/(M[i, pos2idx(j,0)]))
        end
    end
    
    #productivity and initial biomass conc
    Prod_vec = zeros(TL,Nz+1)
    Cinit = params.input_biomass_concentration

    for i in 1:TL
        for j in 0:Ny
            for k in 0:Nz

                    #salinity
                    Sout[i,pos2idx(j,0)] = Sal(Tout[i,pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)
                    #height
                    Ht[i,pos2idx(j,0)] = Hght(Tout[i,pos2idx(j,0)],Vavg[i, pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)
                    #increments in z-direction
                    dz_v[i,pos2idx(j,0)] = dz(Tout[i,pos2idx(j,0)],Vavg[i, pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)

                    #different velocity profiles under laminar and turbulent conditions
                    if Re(Ht[i,pos2idx(j,0)],Tout[i,pos2idx(0,0)],Vavg[i,pos2idx(0,0)]) > 2000
                        V_prof_out[i,pos2idx(j,k)] = Vavg[i,pos2idx(j,0)]
                    else
                        V_prof_out[i,pos2idx(j,k)] = params.velocity_profile_lam(Vavg[i,pos2idx(j,0)],k,Ht[i,pos2idx(j,0)])
                    end

                    Prod_vec[i,pos2idx1(k)] = ((Mout[i,pos2idx(Ny_new,k)]*params.volumetric_flow_rate_o - Cinit*params.volumetric_flow_rate_o)* 24.0* 1000)/(W*L)
                    
     
            end
        end
    end

    #volumetric flowrate and concentration adjustment factor (used to adjust biomass and co2 concentrations due to added water, affects light penetration and CO2 concentration terms in mu)
    Q_out = zeros(TL,Nelements) #volumetric flowrate (m3/hr)
    Q_adj = zeros(TL,Nelements) #concentration adjustment factor

    for i in 1:TL
        for j in 0:Ny
            for k in 0:Nz
                Q_out[i,pos2idx(j,k)] = Q(Tout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i],j)
                Q_adj[i,pos2idx(j,k)] = Q_out[i,pos2idx(0,0)]/Q_out[i,pos2idx(j,0)]
            end
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

    watt_to_umolm2s = 0.425*4.6

    #average light intensity for each wavelength
    I_avg_vec = zeros(301,NelementsT2)

    for l = 1:301
        for i = 1:TL
            for j = 0:Ny
                for k = 0:Nz
                    if k == 0
                        I_avg_vec[l,tpos2idx2(i,j,k)] = GHIout[i]*watt_to_umolm2s*percent_light[l]
                    else
                        I_avg_vec[l,tpos2idx2(i,j,k)] = I_avg_vec[l,tpos2idx2(i,j,k-1)]*exp(-absorbance[l]*Q_adj[i,pos2idx(j,k)]*Mout[i,pos2idx(j,k)]*dz_v[i,pos2idx(j,0)])
                    end
                end
            end
        end
    end

    #average light intensity summed across all wavelengths
    I_avg = zeros(TL,Nelements)

    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
                I_avg[i,pos2idx(j,k)] = sum(I_avg_vec[1:301,tpos2idx2(i,j,k)])*(1/watt_to_umolm2s)
            end
        end
    end

    Y_xphm = 1.38/24 #maximal biomass yield on light, molx/mol
    a_x = 4.81 #absorbance cross section, m2/mol

    #Growth adjustment factor - light
    phiL = zeros(TL,Nelements)
    bm = params.max_biomass_specific_growth_rate
    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
                phiL[i,pos2idx(j,k)] = tanh((Y_xphm*a_x*I_avg[i,pos2idx(j,k)]*1E-06)/(bm*(1/3600)))
            end
        end
    end
   
    #Growth adjustment factor - CO2
    phiCO2 = zeros(TL, Nelements)

    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
                if  CO2_out[i,pos2idx(j,k)]*Q_adj[i,pos2idx(j,k)] <= 0.1782/0.0315 && CO2_out[i,pos2idx(j,k)]*Q_adj[i,pos2idx(j,k)] >= 0.689655
                    phiCO2[i,pos2idx(j,k)] = (0.0261*CO2_out[i,pos2idx(j,k)]*Q_adj[i,pos2idx(j,k)] - 0.018)/0.129383
                elseif CO2_out[i,pos2idx(j,k)]*Q_adj[i,pos2idx(j,k)] > 0.1782/0.0315 && CO2_out[i,pos2idx(j,k)]*Q_adj[i,pos2idx(j,k)] <= 29.67
                    phiCO2[i,pos2idx(j,k)] = (-0.0054*CO2_out[i,pos2idx(j,k)]*Q_adj[i,pos2idx(j,k)] + 0.1602)/0.129383
                else
                    phiCO2[i,pos2idx(j,k)] = 0
                end
            end
        end
    end

    #Adjusted biomass concentration and specific growth rate
    Mout_adj = zeros(TL,Nelements)
    mu_out = zeros(TL,Nz+1)
    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
                mu_out[i,pos2idx1(k)] = params.biomass_specific_growth_rate(Tout[i,pos2idx(Ny_new,k)], Sout[i,pos2idx(Ny_new,0)])*phiL[i,pos2idx(Ny_new,k)]*phiCO2[i,pos2idx(j,k)]-0.003621 
                Mout_adj[i,pos2idx(j,k)] = Mout[i,pos2idx(j,k)]*Q_adj[i,pos2idx(j,k)]
            end
        end
    end

    #Function averages across vertical slice (from just below surface)
    P(C) = Statistics.mean(C[:, pos2idx1(1:Nz)], dims = 2)

    #Averages across last 100 hrs of run:
    Average_Continuous_Productivity = Statistics.mean(P(Prod_vec)[max(1,TL-100):TL]) #g/m2/day
    Average_Spec_Growth = Statistics.mean(P(mu_out)[max(1,TL-100):TL]) #hr-1
    Outlet_Height = Ht[TL,pos2idx(Ny_new,0)] #m

    #DIC and CO2
    DICout_adj = zeros(TL, Nelements)
    CO2out_adj = zeros(TL, Nelements)

    #DIC concentration is adjusted bc of dilution
    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
                DICout_adj[i,pos2idx(j,k)] = DIC_out[i,pos2idx(j,k)]*Q_adj[i,pos2idx(j,k)]
                CO2out_adj[i,pos2idx(j,k)] = CO2_out[i,pos2idx(j,k)]*Q_adj[i,pos2idx(j,k)]
            end
        end
    end

    #prevents DIC_out from exiting interpolation bounds
    DICout_adj = max.(DICout_adj,1E-09)
    DICout_adj = min.(DICout_adj,70000)

    #Average DIC at outlet, avg across last 100 hrs
    Average_DIC = 0
    for j = 1:Nz
       Average_DIC = Average_DIC + Statistics.mean(DICout_adj[max(1,TL-100):TL,pos2idx(Ny_new,j)])
    end
   
    Average_DIC =  Average_DIC/(Nz)

    pH = zeros(TL,Nelements)
    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
                pH[i,pos2idx(j,k)] = params.pH_interp(DICout_adj[i,pos2idx(j,k)]) #takes input in umol/L
            end
        end
    end

    #Average pH at outlet, avg across last 100 hrs
    Avg_pH_out = 0
    for j = 1:Nz
        Avg_pH_out = Avg_pH_out + 10^(Statistics.mean(-pH[max(1,TL-100):TL,pos2idx(Ny_new,j)]))
    end

    Avg_pH_out = -log10(Avg_pH_out/(Nz))

    #Average biomass conc at outlet, avg across last 100 hrs
    Avg_BM_out = 0
    for j = 1:Nz
        Avg_BM_out = Avg_BM_out + Statistics.mean(Mout_adj[max(1,TL-100):TL,pos2idx(Ny_new,j)])
    end

    Avg_BM_out = Avg_BM_out/(Nz)

    return Average_Continuous_Productivity, Average_DIC, Outlet_Height, Avg_pH_out, Avg_BM_out, Average_Spec_Growth, Cum_Evap_Loss
end
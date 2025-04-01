
using CSV, Tables

function Plot_Biomass_Profile(Mout, CO2_out, DIC_out,Tout, A_out, T, params, filesuffix)


    Nx = params.num_odes_x
    Ny = params.num_odes_y
    L = params.reactor_length
    RL = params.real_length #Corresponds to 'actual length' of river reactor (outlet value taken at this position)
    Nx_new = trunc(Int,(RL/L)*Nx) #Index corresponding to 'actual length'
    W = params.reactor_width                    
    TL = length(T) #number of elements in timestep
    dx = L/Nx

    #Number of total elements in vector
    Nelements = (Nx+1)
    NelementsT2 = (Nelements)*(TL+1) #accounts for length, height, time

    #Reindexing of X-Dim position in reactor to be contained in 1D vector
    pos2idx(z) = z.+1 #1-D to 1-D
    pos2idx2(y,z) = (y.+1) .+ z.*(Nx.+1) #2-D to 1-D
    
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
    
    P_a(Tamb) = params.saturated_vapor_pressure_water_air(Tamb) #Pa
    Pa_out = zeros(TL, 1) #Pa
    for i in 1:TL
        Pa_out[i] = P_a(Tambout[i]) #pressure at ambient temp, Pa
    end

    ## Evaporation Loss Calculations

    Evap_C = zeros(TL, Nx+1) #evaporation constant, dimless
    X_air = zeros(TL, Nx+1) #ambient humidity ratio, dimless kg H2O/kg dry air
    P_atm = params.reference_pressure #Pa

    for i in 1:TL
        Evap_C[i,pos2idx(0:Nx)] .= 25 + 19*WNDSPDout[i]
        X_air[i,pos2idx(0:Nx)] .= 0.62198*(Pa_out[i]*(RHout[i]/100))/(P_atm-(Pa_out[i]*(RHout[i]/100)))
    end
   
    X_surface = zeros(TL, Nx+1) #surface humidity ratio, dimless
    M_Evap = zeros(TL, Nx+1) #evaporation mass flux, kg H2O/m2/hr
    for i in 1:TL
        for j in 0:Nx
           
            X_surface[i,pos2idx(j)] = 0.62198*(Pa_out[i]/(P_atm-Pa_out[i])) #kg H2O/kg dry air
            M_Evap[i,pos2idx(j)] = Evap_C[i,pos2idx(j)]*(X_surface[i,pos2idx(j)]-X_air[i, pos2idx(j)]) #kg H2O/m2-hr
        end
    end

    Evap_Loss_T = zeros(TL+1) #loss of water due to evaporation for each time increment
    for i = 1:TL
        Evap_Loss_T[i] = W*RL*M_Evap[i,pos2idx(Nx_new)]*(T[i]-T[max(i-1, 1)]) #selection of j index does not matter, M_evap is constant across reactor
    end

    Cum_Evap_Loss = sum(Evap_Loss_T) #sums across all time increments


    #Functions for height, vertical increment, volumetric flowrate, reynolds number, salinity
    
    Sal(Area) = params.salinity(Area)

    #Initialize vectors
    Sout = zeros(TL, Nx+1)
    Ht = zeros(TL,Nx+1)

    #Vavg determined assuming conservation of kinetic energy (evaporated water has 0 kinetic energy)
    Vavg = zeros(TL, Nelements) #average velocity (m/hr)
    

    #productivity and initial biomass conc
    Prod_vec = zeros(TL,Ny+1)
    Cinit = params.input_biomass_concentration

    alpha = params.alpha
    m = params.m

    #volumetric flowrate and concentration adjustment factor (used to adjust biomass and co2 concentrations due to added water, affects light penetration and CO2 concentration terms in mu)
    Q_out = zeros(TL,Nelements) #volumetric flowrate (m3/hr)

    for i in 1:TL
        for j in 0:Nx

                    Vavg[i,pos2idx(j)] = alpha*A_out[i,pos2idx(j)]^(m-1)
                    Q_out[i,pos2idx(j)] = Vavg[i,pos2idx(j)]*A_out[i,pos2idx(j)]
                    Ht[i,pos2idx(j)] = A_out[i,pos2idx(j)]/W

                    #salinity
                    Sout[i,pos2idx(j)] = Sal(A_out[i,pos2idx(j)])
                    #height
                    Ht[i,pos2idx(j)] = A_out[i,pos2idx(j)]/W

                    Prod_vec[i] = ((Mout[i,pos2idx(Nx_new)]*Q_out[i,pos2idx(Nx_new)] - Cinit*params.volumetric_flow_rate_o)* 24.0* 1000)/(W*L)

                    
        end
    end

    #Similarly, integrate 1/V(x)dx over 0-L to find retention time
    RT_T = zeros(TL, Nelements)
    for i = 1:TL
        for j = 1:Nx
            RT_T[i,pos2idx(j)] = (1/Vavg[i,pos2idx(j)]) *(dx) #selection of j index does not matter, M_evap is constant across reactor
        end
    end
    RT = zeros(TL)
    for i = 1:TL
        for j = 0:Nx
            RT[i] = RT[i] + RT_T[i,pos2idx(j)]
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
            for j = 0:Nx
                I_avg_vec[l,pos2idx2(j,i)] = GHIout[i]*watt_to_umolm2s*percent_light[l]*exp(-absorbance[l]*Mout[i,pos2idx(j)]*Ht[i,pos2idx(j)]/2)
            end
        end
    end

    #average light intensity summed across all wavelengths
    I_avg = zeros(TL,Nelements)

    for i = 1:TL
        for j = 0:Nx
                I_avg[i,pos2idx(j)] = sum(I_avg_vec[1:301,pos2idx2(j,i)])*(1/watt_to_umolm2s)
        end
    end

    Y_xphm = 1.38/24 #maximal biomass yield on light, molx/mol
    a_x = 4.81 #absorbance cross section, m2/mol

    #Growth adjustment factor - light
    phiL = zeros(TL,Nelements)
    bm = params.max_biomass_specific_growth_rate
    for i = 1:TL
        for j = 0:Nx
            phiL[i,pos2idx(j)] = tanh((Y_xphm*a_x*I_avg[i,pos2idx(j)]*1E-06)/(bm*(1/3600)))
        end
    end
   
    #Growth adjustment factor - CO2
    phiCO2 = zeros(TL, Nelements)

    for i = 1:TL
        for j = 0:Nx
                if  CO2_out[i,pos2idx(j)] <= 0.1782/0.0315 && CO2_out[i,pos2idx(j)] >= 0.689655
                    phiCO2[i,pos2idx(j)] = (0.0261*CO2_out[i,pos2idx(j)] - 0.018)/0.129383
                elseif CO2_out[i,pos2idx(j)] > 0.1782/0.0315 && CO2_out[i,pos2idx(j)] <= 29.67
                    phiCO2[i,pos2idx(j)] = (-0.0054*CO2_out[i,pos2idx(j)] + 0.1602)/0.129383
                else
                    phiCO2[i,pos2idx(j)] = 0
                end
        end
    end

    #Adjusted biomass concentration and specific growth rate
    mu = zeros(TL,Nx_new+1)
    mu_out = zeros(TL)
    CO2_star(T) = params.C_CO2_star(T)
    k_co2(WNDSPD,T) = params.k_co2(WNDSPD,T)

    CO2_uptake = zeros(TL, Nelements)
    CO2_flux = zeros(TL, Nelements)

    for i = 1:TL
        for j = 1:Nx
            mu[i,pos2idx(j)] = params.biomass_specific_growth_rate(Tout[i,pos2idx(j)], Sout[i,pos2idx(j)])*phiL[i,pos2idx(j)]*phiCO2[i,pos2idx(j)]-0.003621
            mu_out[i] = mu[i,pos2idx(Nx_new)]

            CO2_uptake[i, pos2idx(j)] = -params.co2_per_biomass * mu[i,pos2idx(j)] * ((Mout[i,pos2idx(j)])*1000)
            CO2_flux[i, pos2idx(j)] = k_co2(WNDSPDout[i],Tout[i,pos2idx(j)])*(CO2_star(Tout[i,pos2idx(j)]) - CO2_out[i,pos2idx(j)])/Ht[i,pos2idx(j)]
        end
    end

    


    #Function averages across vertical slice (from just below surface)

    #Averages across last 100 hrs of run:
    Average_Continuous_Productivity = Statistics.mean(Prod_vec[max(1,TL-100):TL]) #g/m2/day
    Average_CO2_Uptake = Average_Continuous_Productivity*params.co2_per_biomass #g Co2/m2/day
    Average_Spec_Growth = Statistics.mean(mu_out[max(1,TL-100):TL]) #hr-1
    Outlet_Height = Statistics.mean(Ht[max(1,TL-100):TL,pos2idx(Nx_new)]) #m
    Average_RT = Statistics.mean(RT[max(1,TL-100):TL])
 

    #Average DIC at outlet, avg across last 100 hrs
    Average_DIC = Statistics.mean(DIC_out[max(1,TL-100):TL,pos2idx(Nx_new)])
 
    #Average biomass conc at outlet, avg across last 100 hrs
    Avg_BM_out = Statistics.mean(Mout[max(1,TL-100):TL,pos2idx(Nx_new)])

    @show Average_Continuous_Productivity

    Mout2D = zeros(Ny+1, Nx_new+1)
    for i in 0:Nx_new
        for j in 0:Ny
            Mout2D[j+1,i+1] = Statistics.mean(Mout[max(1,TL-100):TL, pos2idx(i)])
        end
    end

    Qout2D = zeros(Ny+1,Nx_new+1)
    for i in 0:Nx_new
        for j in 0:Ny
            Qout2D[j+1,i+1] = Vavg[TL,pos2idx(i)]*A_out[TL,pos2idx(i)]
        end
    end

    Cout2D = zeros(Ny+1, Nx_new+1)
    for i in 0:Nx_new
        for j in 0:Ny
            Cout2D[j+1,i+1] = Statistics.mean(CO2_out[max(1,TL-100):TL, pos2idx(i)]) # don't allow negative values
        end
    end

    Tout2D = zeros(Ny+1, Nx_new+1)
    for i in 0:Nx_new
        for j in 0:Ny
            Tout2D[j+1,i+1] = Statistics.mean(Tout[max(1,TL-100):TL, pos2idx(i)])
        end
    end

    Hout2D = zeros(Ny+1, Nx_new+1)
    for i in 0:Nx_new
        for j in 0:Ny
            Hout2D[j+1,i+1] = Statistics.mean(Ht[max(1,TL-100):TL, pos2idx(i)])
        end
    end

    Lout2D = zeros(Ny+1, Nx_new+1)
    for i in 0:Nx_new
        for j in 0:Ny
            Lout2D[j+1,i+1] = Statistics.mean(I_avg[max(1,TL-100):TL, pos2idx(i)])
        end
    end

    muout2D = zeros(Ny+1, Nx_new+1)
    for i in 0:Nx_new
        for j in 0:Ny
            muout2D[j+1,i+1] = Statistics.mean(mu[max(1,TL-100):TL, pos2idx(i)])
        end
    end

    


    Y = LinRange(0,dx*Nx_new,Nx_new+1)
    Z = LinRange(0,Ht[1,pos2idx(0)],Ny+1)


    q = heatmap(Y, Z, Mout2D,
    yflip=true,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
    title="Algae Cell Density (g/L)",
    size=(800,400)
    )
    savefig(q, "BiomassProfile_AlgaeRiverReactor_$filesuffix.ps")
    png("BiomassProfile_AlgaeRiverReactor_$filesuffix")

    q2 = heatmap(Y,Z, Qout2D,
        yflip=true,
        c=cgrad([:blue, :white,:red, :yellow]),
        xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
        title="Volumetric Flowrate (m^3/hr)",
        size=(800,400)
    )
    savefig(q2, "VolFlrProfile_AlgaeRiverReactor_$filesuffix.ps")
    png("VolFlrProfile_AlgaeRiverReactor_$filesuffix")
    CSV.write("Q_HM_$filesuffix.csv", Tables.table(Cout2D), writeheader = false)
    
    q3 = heatmap(Y, Z, (Cout2D),
            yflip=true,
            c=cgrad([:blue, :white,:red, :yellow]),
            xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
            title="Dissolved CO2 (g/m^3)",
            size=(800,400)
            )
    savefig(q3, "CO2Profile_AlgaeRiverReactor_$filesuffix.ps")
    png("CO2Profile_AlgaeRiverReactor_$filesuffix")

    CSV.write("Temp_HM_$filesuffix.csv", Tables.table(Tout2D), writeheader = false)
    q4 = heatmap(Y, Z, Tout2D,
            yflip=true,
            c=cgrad([:blue, :white,:red, :yellow]),
            xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
            title="Temperature (K)",
            size=(800,400)
            )
    savefig(q4, "TemperatureProfile_AlgaeRiverReactor_$filesuffix.ps")
    png("TemperatureProfile_AlgaeRiverReactor_$filesuffix")

    CSV.write("Hght_HM_$filesuffix.csv", Tables.table(Hout2D), writeheader = false)
    q5 = heatmap(Y, Z, Hout2D,
            yflip=true,
            c=cgrad([:blue, :white,:red, :yellow]),
            xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
            title="Height (m)",
            size=(800,400)
            )
    savefig(q5, "HeightProfile_AlgaeRiverReactor_$filesuffix.ps")
    png("HeightProfile_AlgaeRiverReactor_$filesuffix")

    CSV.write("Light_HM_$filesuffix.csv", Tables.table(Lout2D), writeheader = false)
    q6 = heatmap(Y, Z, Lout2D,
            yflip=true,
            c=cgrad([:blue, :white,:red, :yellow]),
            xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
            title="Light Intensity (W/m^2)",
            size=(800,400)
            )
    savefig(q6, "LightProfile_AlgaeRiverReactor_$filesuffix.ps")
    png("LightProfile_AlgaeRiverReactor_$filesuffix")

    CSV.write("muout_HM_$filesuffix.csv", Tables.table(muout2D), writeheader = false)
    q7 = heatmap(Y, Z, muout2D,
            yflip=true,
            c=cgrad([:blue, :white,:red, :yellow]),
            xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
            title="Specific Growth Rate (1/hr)",
            size=(800,400)
            )
    savefig(q7, "muProfile_AlgaeRiverReactor_$filesuffix.ps")
    png("muProfile_AlgaeRiverReactor_$filesuffix")

    p1 = plot(T[1:TL],Tout[1:TL, pos2idx(Nx_new)], xlabel = "Time [hr]", ylabel = "Temperature at Outlet (K)", title="Outlet Temp", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    png("Tempout_Time_AlgaeRiverReactor_$filesuffix")
    savefig(p1, "Tempout_Time_AlgaeRiverReactor_$filesuffix.ps")


    p4 = plot(T[1:TL],Prod_vec[1:TL], xlabel = "Time [hr]", ylabel = "Abs Column CO2 Conc (g/m3)", title="Abs Column CO2 Conc Change over Length", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    png("Prod_Time_AlgaeRiverReactor_$filesuffix")
    savefig(p4, "Prod_Time_AlgaeRiverReactor_$filesuffix.ps")
    return Average_Continuous_Productivity, Average_DIC, Outlet_Height, Avg_BM_out, Average_Spec_Growth, Cum_Evap_Loss, Average_RT, Average_CO2_Uptake
end
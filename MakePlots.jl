
using CSV, Tables

function Plot_Biomass_Profile(Mout, CO2_out, DIC_out,Tout, T, params, filesuffix)
    DIC_out = max.(DIC_out,1E-09)
    DIC_out = min.(DIC_out,70000)

    Ny = params.num_odes_y
    Nz = params.num_odes_z
    L = params.reactor_length                   # m
    RL = params.real_length #e.g. real length is 1/10 of reactor length
    length_frac = RL/L
    Ny_new = trunc(Int,(length_frac)*Ny)
    Nelements = (Ny+1) * (Nz+1)
    Nelements1 = Ny + 1
    pos2idx1(z) = z.+1
    tpos2idx(t,y) = Nelements1*t + pos2idx1(y)
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    tpos2idx2(t,y,z) = Nelements*t + pos2idx(y,z)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    
    
    TL = length(T)
    NelementsT2 = (Nelements)*(TL+1)
  
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data
    data_begin = params.data_begin
    bm = params.max_biomass_specific_growth_rate
 
    dy = L/Ny
    W = params.reactor_width                    # m
    H = params.reactor_initial_liquid_level     # m
    Cinit = params.input_biomass_concentration  # kg
    
    Y = LinRange(0,dy*Ny_new,Ny_new+1)
    Z = LinRange(0,H,Nz+1)
 
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

    


    Hght(T,S,RH,WNDSPD,P,x) = params.height(T,S,RH,WNDSPD,P,x) #m height(T,S,R_H,W,P,x)
   
    dz(T,S,RH,WNDSPD,P,x) = Hght(T,S,RH,WNDSPD,P,x)/Nz #m
   
    Q(T,WNDSPD,RH,P,x) = params.volumetric_flow_rate(T,WNDSPD,RH,P,x)      # m^3/hour v
    Re(H,T,V) = params.reynolds_number(H,T,V)
  
    Vavg = zeros(TL, Nelements1)
    M = zeros(TL, Nelements1)

    for t = 1:TL
        for i = 0:Ny
            M[t,pos2idx(i,0)] = params.mass_o(Tout[t,pos2idx(i,0)])
        end
    end

    #Re_star(H,T) = params.rough_reynolds_number(H,T) #unitless
    Sal(T,RH, WNDSPD,P,x) = params.salinity(T,RH,WNDSPD,P,x) #kg/m3

    #Initialize vectors
    Sout = zeros(TL, Ny+1) #kg/m3
    Ht = zeros(TL,Ny+1) #m
    dz_v = zeros(TL,Ny+1) #m
    Re_out = zeros(TL, Ny+1) #unitless
    Prod_vec = zeros(TL,Nz+1) #g/m2/day
    Cout_vec = zeros(TL,Nz+1) #g/m3
    Vol_vec = zeros(TL)
    mu_out = zeros(TL,Nz+1) #hr-1
    V_prof_out = zeros(TL,Nelements) #m/hr
    Hght_out = zeros(TL,Ny+1) #m

    M = zeros(TL, Nelements)

    for i in 1:TL
        M[i, pos2idx(0,0)] = W*dy*H*params.density_solution(Tout[i,pos2idx(0,0)], 35)
        Vavg[i, pos2idx(0,0)] = params.volumetric_flow_rate_o/(W*H)
    end

    for i in 1:TL
        for j in 1:Ny
            M[i, pos2idx(j,0)] = M[i, pos2idx(j-1,0)] - dy*(1/Vavg[i, pos2idx(j-1,0)])*params.evaporation_mass_flux(Tout[i,pos2idx(j-1,0)],WNDSPDout[i],RHout[i],Pa_out[i])*dy*W
            Vavg[i, pos2idx(j,0)] = (M[i, pos2idx(j-1,0)]*Vavg[i, pos2idx(j-1,0)])/(M[i, pos2idx(j,0)])
        end
    end

    #Bd = boundary layer height, m
    
    
    for i in 1:TL
        for j in 0:Ny
            for k in 0:Nz


                    Sout[i,pos2idx(j,0)] = Sal(Tout[i,pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)

                    #Vavg[i,pos2idx(j,0)] = params.avg_velocity(Tout[i,pos2idx(j,0)], Sout[i,pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)
                      
                    Ht[i,pos2idx(j,0)] = Hght(Tout[i,pos2idx(j,0)],Vavg[i, pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)
                    dz_v[i,pos2idx(j,0)] = dz(Tout[i,pos2idx(j,0)],Vavg[i, pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)

                    if Re(Ht[i,pos2idx(j,0)],Tout[i,pos2idx(0,0)],Vavg[i,pos2idx(0,0)]) > 2000
                        V_prof_out[i,pos2idx(j,k)] = Vavg[i,pos2idx(j,0)]
                    else
                        V_prof_out[i,pos2idx(j,k)] = params.velocity_profile_lam(Vavg[i,pos2idx(j,0)],k,Ht[i,pos2idx(j,0)])
                    end
                    Re_out[i,pos2idx(j,0)] = Re(Ht[i,pos2idx(j,0)],Tambout[i],Vavg[i,pos2idx1(j)])

                    Prod_vec[i,pos2idx1(k)] = ((Mout[i,pos2idx(Ny_new,k)]*Q(Tout[i,pos2idx(Ny_new,0)],WNDSPDout[i],RHout[i],Pa_out[i],Ny_new) - Cinit*params.volumetric_flow_rate_o)* 24.0* 1000)/(W*L)
                    Cout_vec[i,pos2idx1(k)] = Mout[i,pos2idx(Ny_new,k)]
                    Vol_vec[i] = Q(Tout[i,pos2idx(Ny_new,0)],WNDSPDout[i],RHout[i],Pa_out[i],Ny_new)
                    

                    #((Mout[i,pos2idx(Ny-1,k)]- Cinit)*V_prof_out[i,pos2idx(Ny-1,k)]*(H/Nz)*W) #still uses initial height at beginning because mass is conserved
                    #((Mout[i,pos2idx(Ny-1,k)]* Q(Tout[i,pos2idx(Ny-1,0)],WNDSPDout[i],RHout[i],Pa_out[i],Ny-1)- Cinit* Q(Tout[i,pos2idx(0,0)],WNDSPDout[i],RHout[i],Pa_out[i],0))* 24.0* 1000)/(W*L)
     
            end
        end
    end

    

    #Prod_vec_2 = zeros(TL)
    #for i in 1:TL
        #for k in 0:Nz
            #Prod_vec_2[i] = Prod_vec[i,pos2idx1(k)] + Prod_vec_2[i]
       # end
    #end

    #Prod_out = (Prod_vec_2.* 24.0.* 1000)./(W*L)
    
    Q_out = zeros(TL,Nelements)
    for i in 1:TL
        for j in 0:Ny
            for k in 0:Nz
                Q_out[i,pos2idx(j,k)] = Q(Tout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i],j)
            end
        end
    end


    Hght_out = zeros(TL,Ny+1) #m
    for i in 1:TL
        for j in 0:Ny

            Hght_out[i,pos2idx(j,0)] = Hght(Tout[i,pos2idx(j,0)],Vavg[i, pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)
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

    Conc_New(M_biomass_z,z) = (sum(M_biomass_z[1:z])/z)
    C_new = zeros(TL,Nelements)
    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
                C_new[i,pos2idx(j,k)] = Conc_New(Mout[i,pos2idx(j,0:Nz)],k+1)
            end
        end
    end

    watt_to_umolm2s = 0.425*4.6

    I_avg_vec = zeros(301,NelementsT2)

    for l = 1:301
        for i = 1:TL
            for j = 0:Ny
                for k = 0:Nz
                    if k == 0
                        I_avg_vec[l,tpos2idx2(i,j,k)] = GHIout[i]*watt_to_umolm2s*percent_light[l]
                    else
                        I_avg_vec[l,tpos2idx2(i,j,k)] = I_avg_vec[l,tpos2idx2(i,j,k-1)]*exp(-absorbance[l]*Mout[i,pos2idx(j,k)]*dz_v[i,pos2idx(j,0)])
                    end
                end
            end
        end
    end


    I_avg = zeros(TL,Nelements)
    watt_to_umolm2s = 0.425*4.6

    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
                I_avg[i,pos2idx(j,k)] = sum(I_avg_vec[1:301,tpos2idx2(i,j,k)])*(1/watt_to_umolm2s)
            end
        end
    end

    Y_xphm = 1.38/24 #maximal biomass yield on light, molx/mol
    a_x = 4.81 #absorbance cross section, m2/mol
    phiL = zeros(TL,Nelements)

    phiCO2 = zeros(TL, Nelements)

    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
                if  CO2_out[i, pos2idx(j,k)] <= 0.1782/0.0315
                    phiCO2[i, pos2idx(j,k)] = (0.0261*CO2_out[i, pos2idx(j,k)] - 0.018)/0.129383
                else
                    phiCO2[i, pos2idx(j,k)] = (-0.0054*CO2_out[i, pos2idx(j,k)] + 0.1602)/0.129383
                end
            end
        end
    end


    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
                phiL[i,pos2idx(j,k)] = tanh((Y_xphm*a_x*I_avg[i,pos2idx(j,k)]*1E-06)/(bm*(1/3600))) #umol to mol
                mu_out[i,pos2idx1(k)] = params.biomass_specific_growth_rate(Tout[i,pos2idx(Ny_new,k)], Sout[i,pos2idx(Ny_new,0)], CO2_out[i,pos2idx(Ny_new,k)])*phiL[i,pos2idx(Ny_new,k)]*phiCO2[i, pos2idx(Ny_new,k)]
            end
        end
    end


    (maxval, maxpos) = findmax(Mout[TL,:])
    (Ny_max, Nz_max) = idx2pos(maxpos)

    P(C) = Statistics.mean(C[:, pos2idx1(1:Nz)], dims = 2) #Used for average productivity
    Q_a(C) = Statistics.mean(C[:, pos2idx(0:Ny_new,0)], dims = 2) #Used for average across surface
    R(C) = Statistics.mean(C[:, pos2idx(0:Ny_new,Nz)], dims = 2) #Used for average across floor

    Avg_Conc_y = zeros(Ny_new+1)
    Mout2D = zeros(Nz+1, Ny_new+1)
    for i in 0:Ny_new
        for j in 0:Nz
            Mout2D[j+1,i+1] = Statistics.mean(Mout[max(1,TL-100):TL, pos2idx(i,j)])
        end
    end

    S(C) = Statistics.mean(C[:, pos2idx(Ny_new-1,0:Nz)], dims = 2)
    Avg(C) = Statistics.mean(C[:, 1], dims = 2)
    
    #this will test length = 10 m
    

    p1 = plot(Y,Mout[TL, pos2idx(0:Ny_new,1)] * 1000.0, xlabel = "Length [m]", ylabel = "Algae Cell Density (g/m^3)", title="Algae Surface Cell Density", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p2 = plot(Y,Mout[TL, pos2idx(0:Ny_new,Nz)] * 1000.0, xlabel = "Length [m]", ylabel = "Algae Cell Density (g/m^3)", title="Algae Floor Cell Density", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    #p3 = plot(Z,V_prof_out[TL, pos2idx(Ny_new,0:Nz)]/3600.0, xlabel = "Height [m]", ylabel = "Fluid Velocity (m/s)", title="Velocity Profile at Outlet", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    #p4 = plot(T,Q_a(Re_out), xlabel = "Time [hours]", ylabel = "Reynold's #", title="Average Reynold's # Over Time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p3 = plot(T[max(1,TL-100):TL], P(Cout_vec)[max(1,TL-100):TL], xlabel = "Time [hours]", ylabel = "Algae Cell Density (g/m^3)", title = "Avg Outlet Cell Density", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p4 = plot(T[max(1,TL-100):TL], P(Prod_vec)[max(1,TL-100):TL], xlabel = "Time [hours]", ylabel = "Algae Productivity (g/m^2/day)", title = "Net Continuous Biomass Productivity", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, ylims = (-20,50), grid = false)
    p5 = plot(T[max(1,TL-100):TL], Vol_vec[max(1,TL-100):TL], xlabel = "Time [hours]", ylabel = "Volumetric Flowrate (m3/hr)", title = "Outlet Vol Flr", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    #p7 = plot(T, R(CO2_out), xlabel = "Time [hours]", ylabel = "CO2 (g/m3)", title = "Avg CO2 Conc Over Time (Floor)", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    #p8 = plot(T, GHIout, xlabel = "Time [hours]", ylabel = "Global Horizontal Irradiance (GHI) [W/m^2]", title = "Solar Energy", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p6 =  plot(Y,Avg_Conc_y[1:Ny_new+1] * 1000.0, xlabel = "Length [m]", ylabel = "Average Algae Cell Density (g/m^3)", title="Algae Cell Density Across Length", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p = plot(p1, p2, p3, p4, p5, p6, layout=(6,1), legend=false, size=(1200,1200))
    png("Biomass_AlgaeRiverReactor_$filesuffix")
    savefig(p, "Biomass_AlgaeRiverReactor_$filesuffix.ps")
    @show Average_Continuous_Productivity = Statistics.mean(P(Prod_vec)[max(1,TL-100):TL])
    @show Avg_Cell_Dens = Statistics.mean(P(Cout_vec)[max(1,TL-100):TL])
    @show Average_Out_Vol_flr = Statistics.mean(Vol_vec[max(1,TL-100):TL])

    Average_Spec_Growth = Statistics.mean(P(mu_out)[max(1,TL-100):TL])
    Average_height = Hght_out[TL,pos2idx(Ny_new,0)]
    Average_Dilution = Statistics.mean(R(Q_out)[max(1,TL-100):TL])/(W*L*Average_height)
    Average_RT = L/(Statistics.mean(Q_a(Q_out)[max(1,TL-100):TL])/(Average_height*W)) #hr
    Average_Vol_Flow = Statistics.mean(R(Q_out)[max(1,TL-100):TL])
    Average_Re = Statistics.mean(Avg(Re_out)[max(1,TL-100):TL])


    Topout = zeros(Ny_new+1,1)
    xpos = zeros(Ny_new + 1,1)
    for i in 0:Ny_new
        Topout[i+1,1] = Mout[TL, pos2idx(i,0)]
        xpos[i+1,1] = i*(L/Ny_new)
    end

    Top_table = hcat(Topout, xpos)
    CSV.write("Reactor_Surface.csv", Tables.table(Top_table), writeheader = false)

    Bottomout = zeros(Ny_new+1,1)
    for i in 0:Ny_new
        Bottomout[i+1,1] = Mout[TL, pos2idx(i,Nz)]
    end

    Bottom_table = hcat(Bottomout, xpos)
    CSV.write("Reactor_Floor.csv", Tables.table(Bottom_table), writeheader = false)

    Vel_out = zeros(Nz + 1, 1)
    ypos = zeros(Nz + 1,1)
    for i in 0:Nz
        Vel_out[i+1,1] = V_prof_out[TL, pos2idx(Ny_new,i)]
        ypos[i+1,1] = Hght_out[TL,pos2idx(Ny_new,0)] - (i*(Hght_out[TL,pos2idx(Ny_new,0)]/Nz))
    end

    Vel_table = hcat(Vel_out, ypos)
    CSV.write("Velocity_Out.csv", Tables.table(Vel_table), writeheader = false)

    Light_out = zeros(Nz + 1, 1)
    for i in 0:Nz
        Light_out[i+1,1] = I_avg[TL,pos2idx(Ny_new,i)]
    end

    L_table = hcat(Light_out, ypos)
    CSV.write("Light_Out.csv", Tables.table(L_table), writeheader = false)

    Prod_T = zeros(TL)
    Time = zeros(TL)
    for i in 1:TL
        Prod_T[i] = P(Prod_vec)[i] #avg at outlet
        Time[i] = i
    end

    Prod_TTable = hcat(Prod_T, Time)
    CSV.write("ProdT_Out.csv", Tables.table(Prod_TTable), writeheader = false)

    Temp_T = zeros(TL)
    for i in 1:TL
        Temp_T[i] = S(Tout)[i]
    end

    Temp_TTable = hcat(Temp_T, Time)
    CSV.write("TempT_Out.csv", Tables.table(Temp_TTable), writeheader = false)


    S_T = zeros(TL)
    for i in 1:TL
        S_T[i] = Sout[i, Ny_new] #avg at outlet
    end

    S_TTable = hcat(S_T, Time)
    CSV.write("ST_Out.csv", Tables.table(S_TTable), writeheader = false)




    # Mout2D(z, y), z = 0 is the surface and z = H is the bottom. For plotting, we will invert the z-scale so that z = 0 is the bottom and z = H is the surface.
  
    Lout2D = zeros(Nz+1,Ny_new+1)
    for i in 0:Ny_new
        for j in 0:Nz
            Lout2D[j+1,i+1] = I_avg[TL,pos2idx(i,j)]
        end
    end

    Lout2D2 = zeros(Nz+1,Ny_new+1)
    for i in 0:Ny_new
        for j in 0:Nz
            Lout2D2[j+1,i+1] = Statistics.mean(I_avg[max(1,TL-100):TL, pos2idx(i,j)])
        end
    end

 
    Vout2D = zeros(Nz+1,Ny_new+1)
    for i in 0:Ny_new
        for j in 0:Nz
            Vout2D[j+1,i+1] = V_prof_out[TL,pos2idx(i,j)]
        end
    end

    
    Qout2D = zeros(Nz+1, Ny_new+1)

    for i in 0:Ny_new
        for j in 0:Nz
            Qout2D[j+1,i+1] = Statistics.mean(Q_out[max(1,TL-100):TL,pos2idx(i,j)])
        end
    end


    CSV.write("Biomass_Dens_$filesuffix.csv", Tables.table(Avg_Conc_y),writeheader = false)
    CSV.write("Biomass_HM_$filesuffix.csv", Tables.table(Mout2D), writeheader = false)
    CSV.write("Vel_HM.csv", Tables.table(Vout2D), writeheader = false)
    CSV.write("phiL_HM_$filesuffix.csv", Tables.table(Lout2D), writeheader = false)
    CSV.write("phiL_HM_2$filesuffix.csv", Tables.table(Lout2D2), writeheader = false)
    CSV.write("Vol_HM_$filesuffix.csv", Tables.table(Qout2D), writeheader = false)

    q = heatmap(Y, Z, Mout2D,
            yflip=true,
            c=cgrad([:blue, :white,:red, :yellow]),
            xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
            title="Algae Cell Density (g/L)",
            size=(800,400)
            )
    savefig(q, "BiomassProfile_AlgaeRiverReactor_$filesuffix.ps")
    png("BiomassProfile_AlgaeRiverReactor_$filesuffix")
    q1 = heatmap(Y,Z, Vout2D,
            yflip=true,
            c=cgrad([:blue, :white,:red, :yellow]),
            xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
            title="Velocity Profile (m/hr)",
            size=(800,400)
            )
    savefig(q1, "VelocityProfile_AlgaeRiverReactor_$filesuffix.ps")
    png("VelocityProfile_AlgaeRiverReactor_$filesuffix")

    q2 = heatmap(Y,Z, Lout2D,
    yflip=true,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
    title="Light Intensity (W/m2)",
    size=(800,400)
    )

    savefig(q2, "LightProfile_AlgaeRiverReactor_$filesuffix.ps")
    png("LightProfile_AlgaeRiverReactor_$filesuffix")



    Average_ratio = 0
    ratio = zeros(Nz+1)
    for j = 0:Nz
        ratio[j+1] = (CO2_out[TL, pos2idx(Ny_new-1,j)]*(1/params.molecular_weight_co2))/(DIC_out[TL,pos2idx(Ny_new-1,j)])
        Average_ratio = Average_ratio + ratio[j+1]
    end
    
    Average_ratio =  Average_ratio/(Nz+1)

    Average_DIC = 0
    for j = 0:Nz
       Average_DIC = Average_DIC + DIC_out[TL,pos2idx(Ny_new-1,j)]
    end
   
    Average_DIC =  Average_DIC/(Nz+1)

    Average_Height = Hght_out[TL,pos2idx(Ny_new-1,0)]

    pH = zeros(TL,Nelements)
    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
                pH[i,pos2idx(j,k)] = params.pH_interp(min(DIC_out[i,pos2idx(j,k)],70000)) #takes input in umol/L
            end
        end
    end

    Avg_pH_out = 0
    for j = 0:Nz
        Avg_pH_out = Avg_pH_out + 10^(-pH[TL,pos2idx(Ny_new-1,j)])
    end

    Avg_pH_out = -log10(Avg_pH_out/(Nz+1))

    Avg_BM_out = Avg_Conc_y[Ny_new]

    return Average_Continuous_Productivity, Average_ratio, Average_DIC, Average_Height, Avg_pH_out, Avg_BM_out, Average_Spec_Growth
end

function Plot_Temperature_Profile(Tout, T, params, filesuffix)
    Ny = params.num_odes_y
    Nz = params.num_odes_z
    L = params.reactor_length                   # m
    RL = params.real_length #e.g. real length is 1/10 of reactor length
    length_frac = RL/L
    Ny_new = trunc(Int,(length_frac)*Ny)
    Nelements = (Ny+1) * (Nz+1)
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
  
    dy = L/Ny
    W = params.reactor_width                    # m
    H = params.reactor_initial_liquid_level     # m
    Y = LinRange(0,Ny_new*dy,Ny_new+1)
    Z = LinRange(0,H,Nz+1)
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data
    data_begin = params.data_begin
    TL = length(T)
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
    (maxval, maxpos) = findmax(Tout[TL,:])
    (Ny_max, Nz_max) = idx2pos(maxpos)
    O(T) = Statistics.mean(T[:, pos2idx(Ny_new-1,0:Nz)], dims=2); #avg inlet
    I(T) = Statistics.mean(T[:, pos2idx(1,0:Nz)], dims=2); #avg outlet
    S(T) = Statistics.mean(T[:,pos2idx(0:Ny_new,0)], dims=2); #avg surface
    F(T) = Statistics.mean(T[:,pos2idx(0:Ny_new,Nz)], dims=2); #avg floor
    D(T) = O(T)-I(T)
    p1 = plot(T, S(Tout), xlabel = "Time [hours]", ylabel = "Temperature (K)", title = "Average Surface Temperature", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p2 = plot(T, F(Tout), xlabel = "Time [hours]", ylabel = "Temperature (K)", title = "Average Floor Temperature", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p3 = plot(T, O(Tout), xlabel = "Time [hours]", ylabel = "Temperature (K)", title = "Average Outlet Temperature", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p4 = plot(T, I(Tout), xlabel = "Time [hours]", ylabel = "Temperature (K)", title = "Average Inlet Temperature", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p5 = plot(T, D(Tout), xlabel = "Time [hours]", ylabel = "Temperature (K)", title = "Average Temperature Difference (Outlet - Inlet)", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p6 = plot(T, GHIout, xlabel = "Time [hours]", ylabel = "Global Horizontal Irradiance (GHI) [W/m^2]", title = "Solar Energy", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p7 = plot(T, WNDSPDout, xlabel = "Time [hours]", ylabel = "Wind Speed (WNDSPD) [m/s]", title = "Wind Speed", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p8 = plot(T, RHout, xlabel = "Time [hours]", ylabel = "Relative Humidity (RH) [%]", title = "Relative Humidity", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p9 = plot(T, Tambout, xlabel = "Time [hours]", ylabel = "Ambient Temperature (Tamb) [K]", title = "Ambient Temperature", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, layout=(9,1), legend=false, size=(1600,1600))
    png("Thermal_AlgaeRiverReactor")
    savefig(p, "Thermal_AlgaeRiverReactor.ps")
    # Mout(y,z), z = 0 is the surface and z = H is the bottom. For plotting, we will invert the z-scale so that z = 0 is the bottom and z = H is the surface.
    
    Tout2D = zeros(Nz+1, Ny_new+1)
    for i in 0:Ny_new
        for j in 0:Nz
            Tout2D[j+1,i+1] = Statistics.mean(Tout[max(1,TL-100):TL, pos2idx(i,j)])
        end
    end
    CSV.write("Temp_HM_$filesuffix.csv", Tables.table(Tout2D), writeheader = false)
    q = heatmap(Y, Z, Tout2D,
            yflip=true,
            c=cgrad([:blue, :white,:red, :yellow]),
            xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
            title="Temperature (K)",
            size=(800,400)
            )
    savefig(q, "TemperatureProfile_AlgaeRiverReactor_$filesuffix.ps")
    png("TemperatureProfile_AlgaeRiverReactor_$filesuffix")
end
function Plot_CO2_Profile(CO2_out, DIC_out,Tout, T, params, filesuffix)
    DIC_out = max.(DIC_out,1E-09)
    DIC_out = min.(DIC_out,70000)
    
    Ny = params.num_odes_y
    L = params.reactor_length                   # m
    RL = params.real_length #e.g. real length is 1/10 of reactor length
    length_frac = RL/L
    Ny_new = trunc(Int,(length_frac)*Ny)
    Nz = params.num_odes_z
    Nelements = (Ny+1) * (Nz+1)
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    TL = length(T)
    NelementsT2 = (Nelements)*(TL+1)
    tpos2idx2(t,y,z) = Nelements*t + pos2idx(y,z)

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


    dy = L/Ny
    H = params.reactor_initial_liquid_level     # m
    Y = LinRange(0,dy*Ny_new,Ny_new+1)
    Z = LinRange(0,H,Nz+1)

    T_days = T ./ 24.0
    (maxval, maxpos) = findmax(CO2_out[TL,:])
    (Ny_max, Nz_max) = idx2pos(maxpos)

    P_a(Tamb) = params.saturated_vapor_pressure_water_air(Tamb) #Pa
    Pa_out = zeros(TL, 1) #Pa
    for i in 1:TL
        Pa_out[i] = P_a(Tambout[i]) 
    end

    Hght(T,S,RH,WNDSPD,P,x) = params.height(T,S,RH,WNDSPD,P,x)
    dz(T,S,RH,WNDSPD,P,x) = Hght(T,S,RH,WNDSPD,P,x)/Nz #m

    S(C) = Statistics.mean(C[:, pos2idx(Ny_new-1,0:Nz)], dims = 2)

    p1 = plot(Y,CO2_out[TL, pos2idx(0:Ny_new,0)] , xlabel = "Length [m]", ylabel = "Dissolved CO2 (g/m^3)", title="CO2 Surface Concentration", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p2 = plot(Y,CO2_out[TL, pos2idx(0:Ny_new,Nz)] , xlabel = "Length [m]", ylabel = "Dissolved CO2 (g/m^3)", title="CO2 Floor Concentration", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    
    Time = zeros(TL)
    CO2_T = zeros(TL)
    for i in 1:TL
        CO2_T[i] = S(CO2_out)[i] #avg at outlet
        Time[i] = i
    end

    CO2_TTable = hcat(CO2_T, Time)
    CSV.write("CO2T_Out.csv", Tables.table(CO2_TTable), writeheader = false)


    p = plot(p1, p2, layout=(2,1), legend=false, size=(1200,1200))
    png("CO2_AlgaeRiverReactor")
    savefig(p, "CO2_AlgaeRiverReactor.ps")
    # Cout(y,z), z = 0 is the surface and z = H is the bottom. For plotting, we will invert the z-scale so that z = 0 is the bottom and z = H is the surface.
    Cout2D = zeros(Nz+1, Ny_new+1)
    for i in 0:Ny_new
        for j in 0:Nz
            Cout2D[j+1,i+1] = Statistics.mean(CO2_out[max(1,TL-100):TL, pos2idx(i,j)]) # don't allow negative values
        end
    end
    CSV.write("CO2_HM_$filesuffix.csv", Tables.table(Cout2D), writeheader = false)
    
    q = heatmap(Y, Z, (Cout2D),
            yflip=true,
            c=cgrad([:blue, :white,:red, :yellow]),
            xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
            title="Dissolved CO2 (g/m^3)",
            size=(800,400)
            )
    savefig(q, "CO2Profile_AlgaeRiverReactor_$filesuffix.ps")
    png("CO2Profile_AlgaeRiverReactor_$filesuffix")

    pH = zeros(TL, Nelements)

    PO3 = params.PO3 #mol/kg soln, https://resourcewatch.org/data/explore/f1aa9ec7-c3b6-441c-b395-96fc796b7612?section=Discover&selectedCollection=&zoom=2.422253880286214&lat=51.07099144291875&lng=-85.84319789585153&pitch=0&bearing=0&basemap=dark&labels=light&layers=%255B%257B%2522dataset%2522%253A%2522f1aa9ec7-c3b6-441c-b395-96fc796b7612%2522%252C%2522opacity%2522%253A1%252C%2522layer%2522%253A%25221122cdbf-cb73-467a-bb25-ad86ac491136%2522%257D%255D&aoi=&page=1&sort=most-viewed&sortDirection=-1
    Si = params.Si #mol/kg soln, https://plymsea.ac.uk/id/eprint/1451/1/The_determination_of_silicate_in_sea_water.pdf
    NH4 = params.NH4 #mol/kg soln, http://www.nine-esf.org/files/obergurgl/presentations/Woodward.pdf 
    P_atm = params.P_atm #atm

    T_in = params.input_temperature
    S_in = params.salinity_in

    Sout = zeros(NelementsT2)
    Vavg = zeros(NelementsT2)
    M = zeros(NelementsT2)
    DIC_out2 = zeros(NelementsT2)
    T_out2 = zeros(NelementsT2)
    Sal(T,RH, WNDSPD,P,x) = params.salinity(T,RH,WNDSPD,P,x)
    DIC_in = params.DIC_init
    CO2_in = params.co2_init
    mw_co2 = 0.04401*1000 #g/mol
    
    for i =1:TL
        for j = 0:Ny
            for k = 0:Nz
                T_out2[tpos2idx2(i,j,k)] = Tout[i,pos2idx(j,k)]
            end
        end
    end


    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
                M[tpos2idx2(i,j,k)] = params.mass_o(Tout[i,pos2idx(j,0)])
                Sout[tpos2idx2(i,j,k)] = Sal(Tout[i,pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)
            end
        end
    end

    
    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
                pH[i,pos2idx(j,k)] = params.pH_interp(min(DIC_out[i,pos2idx(j,k)],70000)) #takes input in umol/L
            end
        end
    end


    pHout2D = zeros(Nz+1, Ny_new+1)
    for i in 0:Ny_new
        for j in 0:Nz
            pHout2D[j+1,i+1] = Statistics.mean(pH[max(1,TL-100):TL, pos2idx(i,j)])
        end
    end

    CSV.write("pH_HM_$filesuffix.csv", Tables.table(pHout2D), writeheader = false)

    q1 = heatmap(Y, Z, pHout2D,
        yflip=true,
        c=cgrad([:blue, :white,:red, :yellow]),
        xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
        title="pH",
        size=(800,400)
    )
    savefig(q1, "pHProfile_AlgaeRiverReactor_$filesuffix.ps")
    png("pHProfile_AlgaeRiverReactor_$filesuffix")

end

function Plot_Height_Profile(Tout, T, params, filesuffix)
    Ny = params.num_odes_y
    L = params.reactor_length                   # m
    RL = params.real_length #e.g. real length is 1/10 of reactor length
    length_frac = RL/L
    Ny_new = trunc(Int,(length_frac)*Ny)
    Nz = params.num_odes_z
    Nelements = (Ny+1) * (Nz+1)
    Nelements1 = Ny+1
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    pos2idx1(z) = z.+1
    tpos2idx(t,y) = Nelements1*t + pos2idx1(y)
    
    TL = length(T)

    NelementsT = (Ny+1)*(TL+1)
    

    dy = L/Ny
    W = params.reactor_width                    # m
    H = params.reactor_initial_liquid_level     # m
    Y = LinRange(0,Ny_new*dy,Ny_new+1)
    Z = LinRange(0,H,Nz+1)

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
    
    Evap_C = zeros(TL, Ny+1) #evaporation constant, dimless
    X_air = zeros(TL, Ny+1) #ambient humidity ratio, dimless kg H2O/kg dry air
    X_surface = zeros(TL, Ny+1) #surface humidity ratio, dimless
    M_Evap = zeros(TL, Ny+1) #evaporation mass flux, kg H2O/m2/hr
    T_new = zeros(TL,Ny+1) #time with position index, hr
    Cum_Evap_Loss_vec = zeros(TL, Ny_new+1) #cumulative loss of water due to evaporation, kg
    Cum_Evap_Loss = zeros(TL+1)

    P_atm = params.reference_pressure #Pa
    P_w = params.saturated_vapor_pressure_water_air #Pa

    P_a(Tamb) = params.saturated_vapor_pressure_water_air(Tamb) # Pa

    Pa_out = zeros(TL, 1) #Pa
    for i in 1:TL
        Pa_out[i] = P_a(Tambout[i]) #Pa
    end
    
    for i in 1:TL
        Evap_C[i,pos2idx(0:Ny,0)] .= 25 + 19*WNDSPDout[i]
        X_air[i,pos2idx(0:Ny,0)] .= 0.62198*(P_w(Tambout[i])*(RHout[i]/100))/(P_atm-(P_w(Tambout[i])*(RHout[i]/100)))
    end
   

    for i in 1:TL
        for j in 0:Ny
           
            X_surface[i,pos2idx(j,0)] = 0.62198*(P_w(Tout[i,pos2idx(j,0)])/(P_atm-P_w(Tout[i,pos2idx(j,0)]))) #kg H2O/kg dry air
            M_Evap[i,pos2idx(j,0)] = Evap_C[i,pos2idx(j,0)]*(X_surface[i,pos2idx(j,0)]-X_air[i, pos2idx(j,0)]) #kg H2O/m2-hr
        end
    end

    for i in 1:TL
        for j in 0:Ny_new
            Cum_Evap_Loss_vec[i,pos2idx(j,0)] = W*dy*M_Evap[i,pos2idx(j,0)] #kg H20/hr
        end
    end
    T_new = zeros(TL+1)
    for i = 1:TL
        T_new[i+1] = T[i]
    end

    Cum_Evap_Loss_T = zeros(TL)
    for i = 1:TL
        for j = 0:Ny_new
            Cum_Evap_Loss_T[i] = Cum_Evap_Loss_T[i] + Cum_Evap_Loss_vec[i,pos2idx(j,0)] #kg H20/hr
        end
        Cum_Evap_Loss[i+1] = Cum_Evap_Loss[i] + Cum_Evap_Loss_T[i]*(T_new[i+1]-T_new[i])
    end

    
    popfirst!(Cum_Evap_Loss)

    VolFLR_Lost = zeros(TL, Ny+1)
    density_water(T) = params.density_water(T)

    for i in 1:TL
        for j in 0:Ny
            
            VolFLR_Lost[i,pos2idx(j,0)] = M_Evap[i,pos2idx(j,0)]*(W/density_water(Tout[i,pos2idx(j,0)]))*L
        end
    end
    Vavg = zeros(TL, Nelements1)
    M = zeros(TL, Nelements1)
 
    for t = 1:TL
        for i = 0:Ny
            M[t,pos2idx(i,0)] = params.mass_o(Tout[t,pos2idx(i,0)])
        end
    end

    Q(T,WNDSPD,RH,P,x) = params.volumetric_flow_rate(T,WNDSPD,RH,P,x)      # m^3/hour
    Re(H,T,V) = params.reynolds_number(H,T,V)
    Sal(T,RH, WNDSPD,P,x) = params.salinity(T,RH,WNDSPD,P,x) #kg/kg sw

    Hght(T,S,RH,WNDSPD,P,x) = params.height(T,S,RH,WNDSPD,P,x) #m
    dz(T,S,RH,WNDSPD,P,x) = Hght(T,S,RH,WNDSPD,P,x)/Nz #m

    Ht = zeros(TL,Ny+1) #m
    dz_v = zeros(TL,Ny+1) #m
    Re_out = zeros(TL, Ny+1) #dimless
    Sout = zeros(TL, Ny+1) #kg/m3
    V_prof_out = zeros(TL,Nelements) #m/hr
    #Bd = boundary layer height, m

    M = zeros(TL, Nelements)

    for i in 1:TL
        M[i, pos2idx(0,0)] = W*dy*H*params.density_solution(Tout[i,pos2idx(0,0)], 35)
        Vavg[i, pos2idx(0,0)] = params.volumetric_flow_rate_o/(W*H)
    end

    for i in 1:TL
        for j in 1:Ny
            M[i, pos2idx(j,0)] = M[i, pos2idx(j-1,0)] - dy*(1/Vavg[i, pos2idx(j-1,0)])*params.evaporation_mass_flux(Tout[i,pos2idx(j-1,0)],WNDSPDout[i],RHout[i],Pa_out[i])*dy*W
            Vavg[i, pos2idx(j,0)] = (M[i, pos2idx(j-1,0)]*Vavg[i, pos2idx(j-1,0)])/(M[i, pos2idx(j,0)])
        end
    end
    
    for i in 1:TL
        for j in 0:Ny
            for k in 0:Nz
                    
                    Sout[i,pos2idx(j,0)] = Sal(Tout[i,pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)

                    #Vavg[i,pos2idx(j,0)] = params.avg_velocity(Tout[i,pos2idx(j,0)], Sout[i,pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)
                      
                    Ht[i,pos2idx(j,0)] = Hght(Tout[i,pos2idx(j,0)],Vavg[i, pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)
                   
                    dz_v[i,pos2idx(j,0)] = dz(Tout[i,pos2idx(j,0)],Vavg[i, pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)
                     
                    if Re(Ht[i,pos2idx(j,0)],Tout[i,pos2idx(0,0)],Vavg[i,pos2idx(0,0)]) > 2000
                        V_prof_out[i,pos2idx(j,k)] = Vavg[i,pos2idx(j,0)]
                    else
                        V_prof_out[i,pos2idx(j,k)] = params.velocity_profile_lam(Vavg[i,pos2idx(j,0)],k,Ht[i,pos2idx(j,0)])
                    end
                    
                    Re_out[i,pos2idx(j,0)] = Re(Ht[i,pos2idx(j,0)],Tambout[i],Vavg[i,pos2idx1(j)])

            end
        end
    end


    Hght_out = zeros(TL,Ny+1) #m

    for i in 1:TL
        for j in 0:Ny

            Hght_out[i,pos2idx(j,0)] = Hght(Tout[i,pos2idx(j,0)],Vavg[i, pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)
        end
    end
    
    Hout2D = zeros(Nz+1, Ny_new+1)
    for i in 0:Ny_new
        for j in 0:Nz
            Hout2D[j+1,i+1] = max(0, Hght_out[TL, pos2idx(i,0)]) # don't allow negative values
        end
    end


    CSV.write("Height_HM.csv", Tables.table(Hout2D), writeheader = false)

    q1 = heatmap(Y, Z, Hout2D,
        yflip=true,
        c=cgrad([:blue, :white,:red, :yellow]),
        xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
        title="Height [m]",
        size=(800,400)
    )
    savefig(q1, "heightProfile_AlgaeRiverReactor_$filesuffix.ps")
    png("heightProfile_AlgaeRiverReactor_$filesuffix")
  
    T_days = T ./ 24.0
    (maxval, maxpos) = findmax(Hght_out[TL,:])
    (Ny_max, Nz_max) = idx2pos(maxpos)
    S(Hgt) = Statistics.mean(Hgt[:,pos2idx(0:Ny_new,0)], dims=2) #mean height across reactor, m

  
    p1 = plot(Y,Hght_out[TL, pos2idx(0:Ny_new,0)], xlabel = "Length [m]", ylabel = "Height [m]", title="Height vs Length", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p2 = plot(T,S(Hght_out), xlabel = "Time [hours]", ylabel = "Average Height [m]", title="Average height over time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p3 = plot(T,S(M_Evap), xlabel = "Time [hours]", ylabel = "Average Evaporation Rate [kg/m2-hr]", title="Average evaporation rate over time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p4 = plot(T, Cum_Evap_Loss, xlabel = "Time [hours]", ylabel = "Average Cumulative Water Loss [kg]", title="Average cumulative water loss over time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p = plot(p1, p2, p3, p4, layout=(4,1), legend=false, size=(1200,1200))
    png("Height_AlgaeRiverReactor_$filesuffix")
    savefig(p, "Height_AlgaeRiverReactor_$filesuffix.ps")
    # Cout(y,z), z = 0 is the surface and z = H is the bottom. For plotting, we will invert the z-scale so that z = 0 is the bottom and z = H is the surface.
    
    return Cum_Evap_Loss[TL]
end

function Plot_Salinity_Profile(Tout, T, params, filesuffix)
    Ny = params.num_odes_y
    L = params.reactor_length                   # m
    RL = params.real_length #e.g. real length is 1/10 of reactor length
    length_frac = RL/L
    Ny_new = trunc(Int,(length_frac)*Ny)
    Nz = params.num_odes_z
    Nelements = (Ny+1) * (Nz+1)
    Nelements1 = Ny + 1
    pos2idx1(z) = z.+1
    tpos2idx(t,y) = Nelements1*t + pos2idx1(y)
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    TL = length(T)
    NelementsT = (Ny+1)*(TL+1)

    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data
    data_begin = params.data_begin

    dy = L/Ny
    H = params.reactor_initial_liquid_level     # m

    
    Y = LinRange(0,dy*Ny_new,Ny_new+1)
    Z = LinRange(0,H,Nz+1)
   

    
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
        Pa_out[i] = P_a(Tambout[i])
    end

    Vavg = zeros(TL, Nelements1)
    M = zeros(TL, Nelements1)

    for t = 1:TL
        for i = 0:Ny
            M[t,pos2idx(i,0)] = params.mass_o(Tout[t,pos2idx(i,0)])
        end
    end


    Q(T,WNDSPD,RH,P,x) = params.volumetric_flow_rate(T,WNDSPD,RH,P,x)      # m^3/hour
    Re(H,T,V) = params.reynolds_number(H,T,V)
    Sal(T,RH, WNDSPD,P,x) = params.salinity(T,RH,WNDSPD,P,x) #kg/m3

    Hght(T,S,RH,WNDSPD,P,x) = params.height(T,S,RH,WNDSPD,P,x) #m
    dz(T,S,RH,WNDSPD,P,x) = Hght(T,S,RH,WNDSPD,P,x)/Nz #m

    Ht = zeros(TL,Ny+1) #m
    dz_v = zeros(TL,Ny+1) #m
    Re_out = zeros(TL, Ny+1) #dimless
    Sout = zeros(TL, Ny+1) #kg/m3
    V_prof_out = zeros(TL,Nelements) #m/hr



    for i in 1:TL
        for j in 0:Ny
            for k in 0:Nz
                    
                    Sout[i,pos2idx(j,0)] = Sal(Tout[i,pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)

                    Vavg[i,pos2idx1(j)] = params.avg_velocity(Tout[i,pos2idx(j,0)], Sout[i,pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)
                      
                    Ht[i,pos2idx(j,0)] = Hght(Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)
                   
                    dz_v[i,pos2idx(j,0)] = dz(Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],RHout[i],WNDSPDout[i],Pa_out[i],j)
                     
                    if Re(Ht[i,pos2idx(j,0)],Tout[i,pos2idx(0,0)],Vavg[i,pos2idx(0,0)]) > 2000
                        V_prof_out[i,pos2idx(j,k)] = Vavg[i,pos2idx(j,0)]
                    else
                        V_prof_out[i,pos2idx(j,k)] = params.velocity_profile_lam(Vavg[i,pos2idx(j,0)],k,Ht[i,pos2idx(j,0)])
                    end
            
                    Re_out[i,pos2idx(j,0)] = Re(Ht[i,pos2idx(j,0)],Tambout[i],Vavg[i,pos2idx1(j)])

            end
        end
    end

    Sout2D = zeros(Nz+1, Ny_new+1)
    for i in 0:Ny_new
        for j in 0:Nz
            Sout2D[j+1,i+1] = Statistics.mean(Sout[max(1,TL-100):TL, pos2idx(i,0)])
        end
    end


    CSV.write("Sal_HM.csv", Tables.table(Sout2D), writeheader = false)

    q1 = heatmap(Y, Z, Sout2D,
        yflip=true,
        c=cgrad([:blue, :white,:red, :yellow]),
        xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
        title="Salinity (ppt)",
        size=(800,400)
    )
    savefig(q1, "SalinityHM_AlgaeRiverReactor_$filesuffix.ps")
    png("SalinityHM_AlgaeRiverReactor_$filesuffix")


    

    P(S) = Statistics.mean(S[:,pos2idx(0:Ny_new,0)], dims = 2) #mean salinity across reactor, kg/m3


    p1 = plot(Y,Sout[TL, pos2idx(0:Ny_new,0)], xlabel = "Length [m]", ylabel = "Salinity (ppt)", title="Reactor Salinity with Length", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p2 = plot(T, P(Sout), xlabel = "Time [hours]", ylabel = "Average Salinity (ppt)", title = "Average Reactor Salinity over Time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)

    p = plot(p1, p2, layout=(2,1), legend=false, size=(1200,1200))
    png("Salinity_AlgaeRiverReactor_$filesuffix")
    savefig(p, "Salinity_AlgaeRiverReactor_$filesuffix.ps")
end
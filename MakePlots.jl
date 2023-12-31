
using CSV, Tables
function Plot_Biomass_Profile(Mout, CO2_out, Tout, T, params, filesuffix)
    Ny = params.num_odes_y
    Nz = params.num_odes_z
    Nelements = (Ny+1) * (Nz+1)
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    pos2idx1(z) = z.+1
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    TL = length(T)
    @show TL
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data
    data_begin = params.data_begin

    L = params.reactor_length                   # m
    W = params.reactor_width                    # m
    H = params.reactor_initial_liquid_level     # m
    Cinit = params.input_biomass_concentration  # kg/m^3
    Cmax = params.max_biomass_concentration     # kg/m^3
    Lmax = params.max_biomass_light_saturation
    dy = L/Ny
    
    Y = LinRange(0,L,Ny+1)
    Z = LinRange(0,H,Nz+1)
    Ieff(Cz, z, H, GHI) = GHI .* 0.45 .*(1 - min(sum(Cz[1:z]*(H/Nz))./ Cmax, 1.0))
    phiL(Cz, z,H, GHI) = Ieff(Cz, z, H, GHI) .* exp(1 - Ieff(Cz, z, H, GHI)./Lmax)./ Lmax

   
    
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
    

    H_o = params.reactor_initial_liquid_level #m
    P_a(Tamb) = params.saturated_vapor_pressure_water_air(Tamb) #Pa
    Pa_out = zeros(TL, 1) #Pa
    for i in 1:TL
        Pa_out[i] = P_a(Tambout[i]) #pressure at ambient temp, Pa
    end

    


    Hght(x,T,S,WNDSPD,RH,P) = params.height(x,T,S,WNDSPD,RH,P,H_o) #m
    dz(x,T,S,WNDSPD,RH,P) = Hght(x,T,S,WNDSPD,RH,P)/Nz #m
    mu(T, S, M_biomass_z, GHI, z, dz, C_co2,Hyd) = params.biomass_specific_growth_rate(T, S, M_biomass_z, GHI, z, dz, C_co2,Hyd) #1/hr
    



    Q(H,V) = params.volumetric_flow_rate(H,V)      # m^3/hour
    Re(H,T,V) = params.reynolds_number(H,T,V)
    Vavg_lam(H,T) = params.average_flow_velocity_lam(H,T) #m/hr
    Vavg(H,B) =  params.average_flow_velocity(H,B) #m/hr
    Re_star(H,T) = params.rough_reynolds_number(H,T) #unitless
    Sal(T,WNDSPD,RH,P,x,V) = params.salinity(T,WNDSPD,RH,P,x,V) #kg/m3

    #Initialize vectors
    Sout = zeros(TL, Ny+1) #kg/m3
    Ht = zeros(TL,Ny+1) #m
    dz_v = zeros(TL,Ny+1) #m
    Re_out = zeros(TL, Ny+1) #unitless
    Prod_out = zeros(TL,Nz+1) #g/m2/day
    V_prof_out = zeros(TL,Nelements) #m/hr
    Hght_out = zeros(TL,Ny+1) #m

    #Bd = boundary layer height, m
    

    for i in 1:TL
        for j in 0:Ny
            for k in 0:Nz
                #V_prof_out[i,pos2idx(j,k)] =  params.velocity_profile_lam(dz(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])*k,Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])
                #if Re(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)],Vavg_lam(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])) > 500
                    #if Re_star(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)]) <= 4
                        #global Bd = params.boundary_layer_height_1(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])
                    #elseif Re_star(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)]) <= 11
                        #global Bd = params.boundary_layer_height_2(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])
                    #elseif Re_star(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)]) <= 70
                        #global Bd = params.boundary_layer_height_3(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])
                    #else 
                        #global Bd = params.boundary_layer_height_4
                    #end
                    #for i in 1:TL
                        #for j in 0:Ny
                            #Sout[i,pos2idx(j,0)] = Sal(Tout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i],j,Vavg(Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Bd))
                        #end
                    #end
                    #for i in 1:TL
                        #for j in 0:Ny
                            #Ht[i,pos2idx(j,0)] = Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])
                        #end
                    #end
                    #for i in 1:TL
                        #for j in 0:Ny
                            #dz_v[i,pos2idx(j,0)] = dz(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])
                        #end
                    #end
                    #if k >= 4*Nz/5
                        #V_prof_out[i,pos2idx(j,k)] = params.velocity_profile_nw(dz(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])*k, Ht[i,pos2idx(j,0)], Bd)
                    #else
                        #V_prof_out[i,pos2idx(j,k)] = params.velocity_profile_w(dz(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])*k,Ht[i,pos2idx(j,0)], Bd)
                    #end

                    #for i = 1:TL
                        #for j = 0:Ny
                            #Re_out[i,pos2idx(j,0)] = Re(Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tambout[i],Vavg(Ht[i,pos2idx(j,0)],Bd))
                        #end
                    #end
                    
                    

                #else
                    
                    Sout[i,pos2idx(j,0)] = Sal(Tout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i],j,Vavg_lam(Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,0)]))
                      
                    Ht[i,pos2idx(j,0)] = Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])
                      
                   
                    dz_v[i,pos2idx(j,0)] = dz(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])
                     
                    V_prof_out[i,pos2idx(j,k)] =  params.velocity_profile_lam(dz_v[i,pos2idx(j,0)]*k,Ht[i,pos2idx(j,0)],Tout[i,pos2idx(j,k)])
            
                    Re_out[i,pos2idx(j,0)] = Re(Ht[i,pos2idx(j,0)],Tambout[i],Vavg_lam(Ht[i,pos2idx(j,0)],Tout[i,pos2idx(j,0)]))

                    Prod_out[i,pos2idx1(k)] = ((Mout[i,pos2idx(Ny-1,k)]- Cinit)* Q(Ht[i,pos2idx(j,0)],Vavg_lam(Ht[i,pos2idx(j,0)],Tout[i,pos2idx(0,k)]))* 24.0* 1000)/(W*L)
                   
                    
                #end
            end
        end
    end

    Q_out = zeros(TL,Nelements)
    for i in 1:TL
        for j in 0:Ny
            for k in 0:Nz
                Q_out[i,pos2idx(j,k)] = Q(Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Vavg_lam(Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,0)]))
            end
        end
    end


    Hght_out = zeros(TL,Ny+1) #m
    for i in 1:TL
        for j in 0:Ny

            Hght_out[i,pos2idx(j,0)] = Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])
        end
    end

    molecular_weight_co2 = params.molecular_weight_co2
    co2_to_M = (1/(molecular_weight_co2*1000*1000)) #converts from g/m3 to mol/L: mol*m3/kg-L
    co2_to_uM = co2_to_M*1E6
    Ca_o = params.initial_ca
    Ca = ones(Nelements,1)*Ca_o

    pH_max = 14
    pH = zeros(TL, Nelements)
    pHCO3 = zeros(TL, Nelements)
    pCO3 = zeros(TL, Nelements)
    TC = zeros(TL, Nelements)

    K1(T,S) = params.K1(T,S)
    K2(T,S) = params.K2(T,S)
    Kcal(T,S) = params.Kcal(T,S)
    
    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
           
                pH[i,pos2idx(j,k)] = min(0.5*(-log10(K2(Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)]))-log10(K1(Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)]))
                -log10(max(CO2_out[i,pos2idx(j,k)]*co2_to_uM,1E-10))+log10(1E6)+log10(Kcal(Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)]))-log10(Ca[pos2idx(j,k)])),Float64(pH_max))
                pHCO3[i,pos2idx(j,k)] = -log10(K1(Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)]))-log10(max(CO2_out[i,pos2idx(j,k)]*co2_to_uM,1E-10))+log10(1E6)+pH[i,pos2idx(j,k)]
                pCO3[i,pos2idx(j,k)] = -log10(K2(Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)]))-pHCO3[i,pos2idx(j,k)]+pH[i,pos2idx(j,k)]
                
                #total carbon concentration is original carbon concentration
                #[CO2_old] = [CO2_new] + [HCO3] + [CO3]
                #dC = [CO2_new]-[CO2_old] = -[HCO3_new]-[CO3]
                TC[i,pos2idx(j,k)] = 10^(-pHCO3[i,pos2idx(j,k)]) + 10^(-pCO3[i,pos2idx(j,k)]) + CO2_out[i,pos2idx(j,k)]*co2_to_M
            end
        end
    end


    L_out = zeros(TL,Nelements)
    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
                L_out[i,pos2idx(j,k)] = Ieff(Mout[i,pos2idx(j,0:Nz)], k, Hght_out[i, pos2idx(j, 0)], GHIout[i])
            end
        end
    end

    (maxval, maxpos) = findmax(Mout[TL,:])
    (Ny_max, Nz_max) = idx2pos(maxpos)

    P(C) = Statistics.mean(C[:, pos2idx1(0:Nz)], dims = 2) #Used for average productivity
    Q_a(C) = Statistics.mean(C[:, pos2idx(0:Ny,0)], dims = 2) #Used for average across surface
    R(C) = Statistics.mean(C[:, pos2idx(0:Ny,Nz)], dims = 2) #Used for average across floor
    S(C) = Statistics.mean(C[:, pos2idx(Ny-1,0:Nz)], dims = 2)
    Avg(C) = Statistics.mean(C[:, 1], dims = 2)
    

    p1 = plot(Y,Mout[TL, pos2idx(0:Ny,0)] * 1000.0, xlabel = "Length [m]", ylabel = "Algae Cell Density (g/m^3)", title="Algae Surface Cell Density", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p2 = plot(Y,Mout[TL, pos2idx(0:Ny,Nz)] * 1000.0, xlabel = "Length [m]", ylabel = "Algae Cell Density (g/m^3)", title="Algae Floor Cell Density", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p3 = plot(Z,V_prof_out[TL, pos2idx(Ny,0:Nz)]/3600.0, xlabel = "Height [m]", ylabel = "Fluid Velocity (m/s)", title="Velocity Profile at Outlet", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p4 = plot(T,Q_a(Re_out), xlabel = "Time [hours]", ylabel = "Reynold's #", title="Average Reynold's # Over Time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p5 = plot(T, P(Prod_out), xlabel = "Time [hours]", ylabel = "Algae Productivity (g/m^2/day)", title = "Net Continuous Biomass Productivity", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p6 = plot(T, GHIout, xlabel = "Time [hours]", ylabel = "Global Horizontal Irradiance (GHI) [W/m^2]", title = "Solar Energy", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p = plot(p1, p2, p3, p4, p5, p6, layout=(6,1), legend=false, size=(1200,1200))
    png("Biomass_AlgaeRiverReactor_$filesuffix")
    savefig(p, "Biomass_AlgaeRiverReactor_$filesuffix.ps")
    Average_Continuous_Productivity = Statistics.mean(P(Prod_out)[max(1,TL-10):TL])
    return Average_Continuous_Productivity
    Average_height = Statistics.mean(Q_a(Hght_out))
    Average_Dilution = Statistics.mean(R(Q_out)[max(1,TL-10):TL])/(W*L*Average_height)
    Average_RT = L/(Statistics.mean(Q_a(Q_out)[max(1,TL-10):TL])/(Average_height*W)) #hr
    Average_Vol_Flow = Statistics.mean(R(Q_out)[max(1,TL-10):TL])
    Average_Re = Statistics.mean(Avg(Re_out)[max(1,TL-10):TL])
    @show Average_RT
    @show Average_Re
    @show Average_Dilution
    @show Average_height
    @show Average_Vol_Flow

    Topout = zeros(Ny+1,1)
    xpos = zeros(Ny + 1,1)
    for i in 0:Ny
        Topout[i+1,1] = Mout[TL, pos2idx(i,0)]
        xpos[i+1,1] = i*(L/Ny)
    end

    Top_table = hcat(Topout, xpos)
    CSV.write("Reactor_Surface.csv", Tables.table(Top_table), writeheader = false)

    Bottomout = zeros(Ny+1,1)
    for i in 0:Ny
        Bottomout[i+1,1] = Mout[TL, pos2idx(i,Nz)]
    end

    Bottom_table = hcat(Bottomout, xpos)
    CSV.write("Reactor_Floor.csv", Tables.table(Bottom_table), writeheader = false)

    Vel_out = zeros(Nz + 1, 1)
    ypos = zeros(Nz + 1,1)
    for i in 0:Nz
        Vel_out[i+1,1] = V_prof_out[TL, pos2idx(Ny,i)]
        ypos[i+1,1] = Hght_out[TL,pos2idx(Ny,0)] - (i*(Hght_out[TL,pos2idx(Ny,0)]/Nz))
    end

    Vel_table = hcat(Vel_out, ypos)
    CSV.write("Velocity_Out.csv", Tables.table(Vel_table), writeheader = false)

    Light_out = zeros(Nz + 1, 1)
    for i in 0:Nz
        Light_out[i+1,1] = L_out[TL,pos2idx(Ny,i)]
    end

    L_table = hcat(Light_out, ypos)
    CSV.write("Light_Out.csv", Tables.table(L_table), writeheader = false)

    Prod_T = zeros(TL)
    Time = zeros(TL)
    for i in 1:TL
        Prod_T[i] = P(Prod_out)[i] #avg at outlet
        Time[i] = i
    end

    Prod_TTable = hcat(Prod_T, Time)
    CSV.write("ProdT_Out.csv", Tables.table(Prod_TTable), writeheader = false)

    Lght_T = zeros(TL)
    for i in 1:TL
        Lght_T[i] = S(L_out)[i] #avg at outlet
    end

    Temp_T = zeros(TL)
    for i in 1:TL
        Temp_T[i] = S(Tout)[i]
    end

    Lght_TTable = hcat(Lght_T, Time)
    CSV.write("LghtT_Out.csv", Tables.table(Lght_TTable), writeheader = false)

    Temp_TTable = hcat(Temp_T, Time)
    CSV.write("TempT_Out.csv", Tables.table(Temp_TTable), writeheader = false)


    S_T = zeros(TL)
    for i in 1:TL
        S_T[i] = Sout[i, Ny] #avg at outlet
    end

    S_TTable = hcat(S_T, Time)
    CSV.write("ST_Out.csv", Tables.table(S_TTable), writeheader = false)

    pH_T = zeros(TL)
    for i in 1:TL
        pH_T[i] = S(pH)[i] #avg at outlet
    end

    pH_TTable = hcat(pH_T, Time)
    CSV.write("pHT_Out.csv", Tables.table(pH_TTable), writeheader = false)



    # Mout2D(z, y), z = 0 is the surface and z = H is the bottom. For plotting, we will invert the z-scale so that z = 0 is the bottom and z = H is the surface.
    Mout2D = zeros(Nz+1, Ny+1)
    for i in 0:Ny
        for j in 0:Nz
            Mout2D[j+1,i+1] = Mout[TL, pos2idx(i,j)]
        end
    end

    Lout2D = zeros(Nz+1,Ny+1)
    for i in 0:Ny
        for j in 0:Nz
            Lout2D[j+1,i+1] = L_out[TL,pos2idx(i,j)]
        end
    end
 
    Vout2D = zeros(Nz+1,Ny+1)
    for i in 0:Ny
        for j in 0:Nz
            Vout2D[j+1,i+1] = V_prof_out[TL,pos2idx(i,j)]
        end
    end
    CSV.write("Biomass_HM.csv", Tables.table(Mout2D), writeheader = false)
    CSV.write("Vel_HM.csv", Tables.table(Vout2D), writeheader = false)
    CSV.write("phiL_HM.csv", Tables.table(Lout2D), writeheader = false)

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
end

function Plot_Temperature_Profile(Tout, T, params, filesuffix)
    Ny = params.num_odes_y
    Nz = params.num_odes_z
    Nelements = (Ny+1) * (Nz+1)
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    L = params.reactor_length                   # m
    W = params.reactor_width                    # m
    H = params.reactor_initial_liquid_level     # m
    Y = LinRange(0,L,Ny+1)
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
    O(T) = Statistics.mean(T[:, pos2idx(Ny-1,0:Nz)], dims=2); #avg inlet
    I(T) = Statistics.mean(T[:, pos2idx(1,0:Nz)], dims=2); #avg outlet
    S(T) = Statistics.mean(T[:,pos2idx(0:Ny,0)], dims=2); #avg surface
    F(T) = Statistics.mean(T[:,pos2idx(0:Ny,Nz)], dims=2); #avg floor
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
    png("Thermal_AlgaeRiverReactor_$filesuffix")
    savefig(p, "Thermal_AlgaeRiverReactor_$filesuffix.ps")
    # Mout(y,z), z = 0 is the surface and z = H is the bottom. For plotting, we will invert the z-scale so that z = 0 is the bottom and z = H is the surface.
    Tout2D = zeros(Nz+1, Ny+1)
    for i in 0:Ny
        for j in 0:Nz
            Tout2D[j+1,i+1] = Tout[TL, pos2idx(i,j)]
        end
    end
    CSV.write("Temp_HM.csv", Tables.table(Tout2D), writeheader = false)
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
function Plot_CO2_Profile(CO2_out, Tout, T, params, filesuffix)
    Ny = params.num_odes_y
    Nz = params.num_odes_z
    Nelements = (Ny+1) * (Nz+1)
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    TL = length(T)

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

    L = params.reactor_length                   # m
    W = params.reactor_width                    # m
    H = params.reactor_initial_liquid_level     # m
    Cinit = params.input_biomass_concentration  # kg/m^3
    Cmax = params.max_biomass_concentration     # kg/m^3
    Y = LinRange(0,L,Ny+1)
    Z = LinRange(0,H,Nz+1)
    dy = L/Ny

    Ieff(Cz, z, Hgt) = (1 - min(sum(Cz[1:z]*(H/Nz))./ Cmax, 1.0))
    phiL(Cz, zlist,Hgt) = [Ieff(Cz, z, Hgt) * exp(1 - Ieff(Cz, z, Hgt)) for z in zlist]
    T_days = T ./ 24.0
    (maxval, maxpos) = findmax(CO2_out[TL,:])
    (Ny_max, Nz_max) = idx2pos(maxpos)



    H_o = params.reactor_initial_liquid_level

    P_a(Tamb) = params.saturated_vapor_pressure_water_air(Tamb) #Pa
    Pa_out = zeros(TL, 1) #Pa
    for i in 1:TL
        Pa_out[i] = P_a(Tambout[i]) 
    end

    Hght(x,T,S,WNDSPD,RH,P) = params.height(x,T,S,WNDSPD,RH,P,H_o) #m
    dz(x,T,S,WNDSPD,RH,P) = Hght(x,T,S,WNDSPD,RH,P)/Nz #m

    S(C) = Statistics.mean(C[:, pos2idx(Ny-1,0:Nz)], dims = 2)

    p1 = plot(Y,CO2_out[TL, pos2idx(0:Ny,0)] , xlabel = "Length [m]", ylabel = "Dissolved CO2 (g/m^3)", title="CO2 Surface Concentration", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p2 = plot(Y,CO2_out[TL, pos2idx(0:Ny,Nz)] , xlabel = "Length [m]", ylabel = "Dissolved CO2 (g/m^3)", title="CO2 Floor Concentration", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    @show 
    
    Time = zeros(TL)
    CO2_T = zeros(TL)
    for i in 1:TL
        CO2_T[i] = S(CO2_out)[i] #avg at outlet
        Time[i] = i
    end

    CO2_TTable = hcat(CO2_T, Time)
    CSV.write("CO2T_Out.csv", Tables.table(CO2_TTable), writeheader = false)


    p = plot(p1, p2, layout=(2,1), legend=false, size=(1200,1200))
    png("CO2_AlgaeRiverReactor_$filesuffix")
    savefig(p, "CO2_AlgaeRiverReactor_$filesuffix.ps")
    # Cout(y,z), z = 0 is the surface and z = H is the bottom. For plotting, we will invert the z-scale so that z = 0 is the bottom and z = H is the surface.
    Cout2D = zeros(Nz+1, Ny+1)
    for i in 0:Ny
        for j in 0:Nz
            Cout2D[j+1,i+1] = max(0, CO2_out[TL, pos2idx(i,j)]) # don't allow negative values
        end
    end
    CSV.write("CO2_HM.csv", Tables.table(Cout2D), writeheader = false)
    
    q = heatmap(Y, Z, (Cout2D),
            yflip=true,
            c=cgrad([:blue, :white,:red, :yellow]),
            xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
            title="Dissolved CO2 (g/m^3)",
            size=(800,400)
            )
    savefig(q, "CO2Profile_AlgaeRiverReactor_$filesuffix.ps")
    png("CO2Profile_AlgaeRiverReactor_$filesuffix")

    molecular_weight_co2 = params.molecular_weight_co2
    co2_to_M = (1/(molecular_weight_co2*1000*1000)) #converts from g/m3 to mol/L: mol*m3/kg-L
    co2_to_uM = co2_to_M*1E6
    Ca_o = params.initial_ca
    Ca = ones(Nelements,1)*Ca_o

    pH_max = 14
    pH = zeros(TL, Nelements)
    pHCO3 = zeros(TL, Nelements)
    pCO3 = zeros(TL, Nelements)
    TC = zeros(TL, Nelements)

    K1(T,S) = params.K1(T,S)
    K2(T,S) = params.K2(T,S)
    Kcal(T,S) = params.Kcal(T,S)

    Sout = ones(TL,Nelements)*26.44
    
    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
           
                pH[i,pos2idx(j,k)] = min(0.5*(-log10(K2(Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)]))-log10(K1(Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)]))
                -log10(max(CO2_out[i,pos2idx(j,k)]*co2_to_uM,1E-10))+log10(1E6)+log10(Kcal(Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)]))-log10(Ca[pos2idx(j,k)])),Float64(pH_max))
                #min(0.5*(-log10(K2(Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)]))-log10(K1(Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)]))-log10(max(CO2_out[i,pos2idx(j,k)]*co2_to_M,1E-10))+log10(Kcal(Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)]))-log10(Ca[pos2idx(j,k)])),Float64(pH_max))
                #pOH[pos2idx(i,j)] = -log10(KW) - pH[i,pos2idx(j,k)]
                pHCO3[i,pos2idx(j,k)] = -log10(K1(Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)]))-log10(max(CO2_out[i,pos2idx(j,k)]*co2_to_uM,1E-10))+log10(1E6)+pH[i,pos2idx(j,k)]
                pCO3[i,pos2idx(j,k)] = -log10(K2(Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)]))-pHCO3[i,pos2idx(j,k)]+pH[i,pos2idx(j,k)]
                
                #total carbon concentration is original carbon concentration
                #[CO2_old] = [CO2_new] + [HCO3] + [CO3]
                #dC = [CO2_new]-[CO2_old] = -[HCO3_new]-[CO3]
                TC[i,pos2idx(j,k)] = 10^(-pHCO3[i,pos2idx(j,k)]) + 10^(-pCO3[i,pos2idx(j,k)]) + CO2_out[i,pos2idx(j,k)]*co2_to_M
            end
        end
    end

    pHout2D = zeros(Nz+1, Ny+1)
    for i in 0:Ny
        for j in 0:Nz
            pHout2D[j+1,i+1] = max(0, pH[TL, pos2idx(i,j)]) # don't allow negative values
        end
    end

    CSV.write("pH_HM.csv", Tables.table(pHout2D), writeheader = false)

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
    Nz = params.num_odes_z
    Nelements = (Ny+1) * (Nz+1)
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    TL = length(T)
    L = params.reactor_length                   # m
    W = params.reactor_width                    # m
    H = params.reactor_initial_liquid_level     # m
    Cinit = params.input_biomass_concentration  # kg/m^3
    Cmax = params.max_biomass_concentration     # kg/m^3
    Y = LinRange(0,L,Ny+1)
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
    Cum_Evap_Loss = zeros(TL, Ny+1) #cumulative loss of water due to evaporation, kg

    P_atm = params.reference_pressure #Pa
    P_w = params.saturated_vapor_pressure_water_air #Pa

    H_o = params.reactor_initial_liquid_level #m

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
        T_new[i,pos2idx(0:Ny,0)] .= T[i]
    end

    for i in 1:TL
        for j in 0:Ny
           
            X_surface[i,pos2idx(j,0)] = 0.62198*(P_w(Tout[i,pos2idx(j,0)])/(P_atm-P_w(Tout[i,pos2idx(j,0)]))) #kg H2O/kg dry air
            M_Evap[i,pos2idx(j,0)] = Evap_C[i,pos2idx(j,0)]*(X_surface[i,pos2idx(j,0)]-X_air[i, pos2idx(j,0)])
            Cum_Evap_Loss[i,pos2idx(j,0)] = W*L*sum(M_Evap[1:i,pos2idx(j,0)]) #kg
        end
    end

    Q(H,V) = params.volumetric_flow_rate(H,V)      # m^3/hour
    Re(H,T,V) = params.reynolds_number(H,T,V) #dimless
    Vavg_lam(H,T) = params.average_flow_velocity_lam(H,T) #m/hr
    Vavg(H,B) =  params.average_flow_velocity(H,B) #m/hr
    Re_star(H,T) = params.rough_reynolds_number(H,T) #dimless
    Sal(T,WNDSPD,RH,P,x,V) = params.salinity(T,WNDSPD,RH,P,x,V) #kg/m3
    Hght(x,T,S,WNDSPD,RH,P) = params.height(x,T,S,WNDSPD,RH,P,H_o) #m
    dz(x,T,S,WNDSPD,RH,P) = Hght(x,T,S,WNDSPD,RH,P)/Nz #m

    Ht = zeros(TL,Ny+1) #m
    dz_v = zeros(TL,Ny+1) #m
    Re_out = zeros(TL, Ny+1) #dimless
    Sout = zeros(TL, Ny+1) #kg/m3
    V_prof_out = zeros(TL,Nelements) #m/hr
    #Bd = boundary layer height, m
    
    for i in 1:TL
        for j in 0:Ny
            for k in 0:Nz
                #V_prof_out[i,pos2idx(j,k)] =  params.velocity_profile_lam(dz(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])*k,Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])
                #if Re(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)],Vavg_lam(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])) > 500
                    #if Re_star(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)]) <= 4
                        #global Bd = params.boundary_layer_height_1(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])
                    #elseif Re_star(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)]) <= 11
                        #global Bd = params.boundary_layer_height_2(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])
                    #elseif Re_star(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)]) <= 70
                        #global Bd = params.boundary_layer_height_3(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])
                    #else 
                        #global Bd = params.boundary_layer_height_4
                    #end
                    #for i in 1:TL
                        #for j in 0:Ny
                            #Sout[i,pos2idx(j,0)] = Sal(Tout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i],j,Vavg(Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Bd))
                        #end
                    #end
                    #for i in 1:TL
                        #for j in 0:Ny
                            #Ht[i,pos2idx(j,0)] = Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])
                        #end
                    #end
                    #for i in 1:TL
                        #for j in 0:Ny
                            #dz_v[i,pos2idx(j,0)] = dz(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])
                        #end
                    #end
                    #if k >= 4*Nz/5
                        #V_prof_out[i,pos2idx(j,k)] = params.velocity_profile_nw(dz(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])*k, Ht[i,pos2idx(j,0)], Bd)
                    #else
                        #V_prof_out[i,pos2idx(j,k)] = params.velocity_profile_w(dz(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])*k,Ht[i,pos2idx(j,0)], Bd)
                    #end

                    #for i = 1:TL
                        #for j = 0:Ny
                            #Re_out[i,pos2idx(j,0)] = Re(Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tambout[i],Vavg(Ht[i,pos2idx(j,0)],Bd))
                        #end
                    #end
                    
                    

                #else
                    
                    Sout[i,pos2idx(j,0)] = Sal(Tout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i],j,Vavg_lam(Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,0)]))
                      
                    Ht[i,pos2idx(j,0)] = Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])
                      
                   
                    dz_v[i,pos2idx(j,0)] = dz(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])
                     
                    V_prof_out[i,pos2idx(j,k)] =  params.velocity_profile_lam(dz_v[i,pos2idx(j,0)]*k,Ht[i,pos2idx(j,0)],Tout[i,pos2idx(j,k)])
            
                    Re_out[i,pos2idx(j,0)] = Re(Ht[i,pos2idx(j,0)],Tambout[i],Vavg_lam(Ht[i,pos2idx(j,0)],Tout[i,pos2idx(j,0)]))
                    
                   
                    
                #end
            end
        end
    end

    K = params.dVavgdx #hr-1
    rho_sol = params.density_solution #kg/m3
    Kterm_out = zeros(TL, Ny+1) #hr-1
    for i in 1:TL
        for j in 0:Ny
            Kterm_out[i,pos2idx(j,0)] = K(Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],RHout[i],WNDSPDout[i],P_w(Tambout[i]))*rho_sol(Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)])
        end
    end


    Hght_out = zeros(TL,Ny+1) #m

    for i in 1:TL
        for j in 0:Ny

            Hght_out[i,pos2idx(j,0)] = Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])
        end
    end
    
    Hout2D = zeros(Nz+1, Ny+1)
    for i in 0:Ny
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
    S(Hgt) = Statistics.mean(Hgt[:,pos2idx(0:Ny,0)], dims=2) #mean height across reactor, m

  
    p1 = plot(Y,Hght_out[TL, pos2idx(0:Ny,0)], xlabel = "Length [m]", ylabel = "Height [m]", title="Height vs Length", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p2 = plot(T,S(Hght_out), xlabel = "Time [hours]", ylabel = "Average Height [m]", title="Average height over time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p3 = plot(T,S(M_Evap), xlabel = "Time [hours]", ylabel = "Average Evaporation Rate [kg/m2-hr]", title="Average evaporation rate over time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p4 = plot(T,S(Kterm_out), xlabel = "Time [hours]", ylabel = "Average Counter-Evaporation Rate [kg/m2-hr]", title="Average counter-evaporation rate over time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p5 = plot(T,S(Cum_Evap_Loss), xlabel = "Time [hours]", ylabel = "Average Cumulative Water Loss [kg]", title="Average cumulative water loss over time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p = plot(p1, p2, p3, p4, p5, layout=(5,1), legend=false, size=(1200,1200))
    png("Height_AlgaeRiverReactor_$filesuffix")
    savefig(p, "Height_AlgaeRiverReactor_$filesuffix.ps")
    # Cout(y,z), z = 0 is the surface and z = H is the bottom. For plotting, we will invert the z-scale so that z = 0 is the bottom and z = H is the surface.
    
end

function Plot_Salinity_Profile(Tout, T, params, filesuffix)
    Ny = params.num_odes_y
    Nz = params.num_odes_z
    Nelements = (Ny+1) * (Nz+1)
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    TL = length(T)
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data
    data_begin = params.data_begin

    L = params.reactor_length                   # m
    W = params.reactor_width                    # m
    H = params.reactor_initial_liquid_level     # m
    Cinit = params.input_biomass_concentration  # kg/m^3
    Cmax = params.max_biomass_concentration     # kg/m^3
    dy = L/Ny
    
    Y = LinRange(0,L,Ny+1)
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
    

    H_o = params.reactor_initial_liquid_level #m
    P_a(Tamb) = params.saturated_vapor_pressure_water_air(Tamb) #Pa
    Pa_out = zeros(TL, 1) #Pa
    for i in 1:TL
        Pa_out[i] = P_a(Tambout[i])
    end


    Q(H,V) = params.volumetric_flow_rate(H,V)       # m^3/hour
    Re(H,T,V) = params.reynolds_number(H,T,V) #dimless
    Vavg_lam(H,T) = params.average_flow_velocity_lam(H,T) #m/hr
    Vavg(H,B) =  params.average_flow_velocity(H,B) #m/hr
    Re_star(H,T) = params.rough_reynolds_number(H,T) #dimless
    Sal(T,WNDSPD,RH,P,x,V) = params.salinity(T,WNDSPD,RH,P,x,V) #kg/m3
    Hght(x,T,S,WNDSPD,RH,P) = params.height(x,T,S,WNDSPD,RH,P,H_o) #m
    dz(x,T,S,WNDSPD,RH,P) = Hght(x,T,S,WNDSPD,RH,P)/Nz #m

    Ht = zeros(TL,Ny+1) #m
    dz_v = zeros(TL,Ny+1) #m
    Re_out = zeros(TL, Ny+1) #dimless
    Sout = zeros(TL, Ny+1) #kg/m3
    V_prof_out = zeros(TL,Nelements) #m/hr

    for i in 1:TL
        for j in 0:Ny
            for k in 0:Nz
                #V_prof_out[i,pos2idx(j,k)] =  params.velocity_profile_lam(dz(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])*k,Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])
                #if Re(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)],Vavg_lam(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])) > 500
                    #if Re_star(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)]) <= 4
                        #global Bd = params.boundary_layer_height_1(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])
                    #elseif Re_star(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)]) <= 11
                        #global Bd = params.boundary_layer_height_2(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])
                    #elseif Re_star(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)]) <= 70
                        #global Bd = params.boundary_layer_height_3(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])
                    #else 
                        #global Bd = params.boundary_layer_height_4
                    #end
                    #for i in 1:TL
                        #for j in 0:Ny
                            #Sout[i,pos2idx(j,0)] = Sal(Tout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i],j,Vavg(Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Bd))
                        #end
                    #end
                    #for i in 1:TL
                        #for j in 0:Ny
                            #Ht[i,pos2idx(j,0)] = Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])
                        #end
                    #end
                    #for i in 1:TL
                        #for j in 0:Ny
                            #dz_v[i,pos2idx(j,0)] = dz(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])
                        #end
                    #end
                    #if k >= 4*Nz/5
                        #V_prof_out[i,pos2idx(j,k)] = params.velocity_profile_nw(dz(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])*k, Ht[i,pos2idx(j,0)], Bd)
                    #else
                        #V_prof_out[i,pos2idx(j,k)] = params.velocity_profile_w(dz(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])*k,Ht[i,pos2idx(j,0)], Bd)
                    #end

                    #for i = 1:TL
                        #for j = 0:Ny
                            #Re_out[i,pos2idx(j,0)] = Re(Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tambout[i],Vavg(Ht[i,pos2idx(j,0)],Bd))
                        #end
                    #end
                    
                    

                #else
                    
                    Sout[i,pos2idx(j,0)] = Sal(Tout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i],j,Vavg_lam(Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,0)]))
                      
                    Ht[i,pos2idx(j,0)] = Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])
                      
                   
                    dz_v[i,pos2idx(j,0)] = dz(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])
                     
                    V_prof_out[i,pos2idx(j,k)] =  params.velocity_profile_lam(dz_v[i,pos2idx(j,0)]*k,Ht[i,pos2idx(j,0)],Tout[i,pos2idx(j,k)])
            
                    Re_out[i,pos2idx(j,0)] = Re(Ht[i,pos2idx(j,0)],Tambout[i],Vavg_lam(Ht[i,pos2idx(j,0)],Tout[i,pos2idx(j,0)]))
                    
                   
                    
                #end
            end
        end
    end

    Sout2D = zeros(Nz+1, Ny+1)
    for i in 0:Ny
        for j in 0:Nz
            Sout2D[j+1,i+1] = max(0, Sout[TL, pos2idx(i,0)]) # don't allow negative values
        end
    end


    CSV.write("Sal_HM.csv", Tables.table(Sout2D), writeheader = false)

    q1 = heatmap(Y, Z, Sout2D,
        yflip=true,
        c=cgrad([:blue, :white,:red, :yellow]),
        xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
        title="Salinity (kg/m3)",
        size=(800,400)
    )
    savefig(q1, "SalinityHM_AlgaeRiverReactor_$filesuffix.ps")
    png("SalinityHM_AlgaeRiverReactor_$filesuffix")
    

    P(S) = Statistics.mean(S[:,pos2idx(0:Ny,0)], dims = 2) #mean salinity across reactor, kg/m3


    p1 = plot(Y,Sout[TL, pos2idx(0:Ny,0)], xlabel = "Length [m]", ylabel = "Salinity (kg/m^3)", title="Reactor Salinity with Length", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p2 = plot(T, P(Sout), xlabel = "Time [hours]", ylabel = "Average Salinity (kg/m^3)", title = "Average Reactor Salinity over Time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)

    p = plot(p1, p2, layout=(2,1), legend=false, size=(1200,1200))
    png("Salinity_AlgaeRiverReactor_$filesuffix")
    savefig(p, "Salinity_AlgaeRiverReactor_$filesuffix.ps")
end
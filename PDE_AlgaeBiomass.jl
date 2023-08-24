cd("C:/Users/slant/OneDrive/Desktop/Julia_2")
function PDE_AlgaeBiomass!(dX, M, Temperature,Height, params, t)

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
    mu_max = params.max_biomass_specific_growth_rate
    Q = params.volumetric_flow_rate
    V_profile(z, H) = params.velocity_profile(z, H, Tamb)
    Dy = params.biomass_diffusion_coefficient_y
    Dz = params.biomass_diffusion_coefficient_z
    L = params.reactor_length
    W = params.reactor_width
    Ny = params.num_odes_y
    Nz = params.num_odes_z

    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]

    dy = L / Ny

    Nelements = (Ny+1) * (Nz+1)
    M = max.(M, 0.0) #don't allow negative values
    dM = zeros( Nelements, 1)
    H_o = params.reactor_initial_liquid_level
    P_a = params.saturated_vapor_pressure_water_air(Tamb) 
    K(T,S) = params.dVavgdx(T,S,RH,WNDSPD,P_a)
    M_Evap(T) = params.evaporation_mass_flux(T,WNDSPD,RH,P_a)
    rho_solution(S) = params.density_solution(Tamb,S)

    Hght(x,T,S) = params.height(x,T,S,WNDSPD,RH,P_a,H_o)
    dz(x,T,S) = Hght(x,T,S)/Nz



    ##convert to mass balance
    Vol(dz_y) = dy*dz_y*W
    Iave(Mz, dz_y, z) = GHI .* 0.45 .* (1 - min(sum(Mz[1:z]/(Vol(dz_y))).*dz_y./Cmax,1.0) ) #
    phiL(Mz, dz_y, z) = Iave(Mz, dz_y, z) .* exp(1 - Iave(Mz, dz_y, z)./Isat) ./ Isat
    mu(Mz, dz_y, z) = mu_max .* phiL(Mz, dz_y, z)

    
        

    # ==========  ALGAE GROWTH ========

    #BC1: C(y,z) at y = 0 is constant [Cin]
    #     dC(y,z) at y = 0 is zero.

    #y = 0, z = 0 ... Nz
    dM[pos2idx(0,0:Nz)] .= 0.0
    #@show size(dC)
    #BC2: molar flux is zero at z = 0, H
    # C[pos2idx(i,-1)] = C[pos2idx(i,0)]
    # C[pos2idx(i,Nz+1)] = C[pos2idx(i,Nz)]
    # C[pos2idx(0,-1)] = C[pos2idx(i,0)]
    # C[pos2idx(0,Nz+1)] = C[pos2idx(i,Nz)]
    # C[pos2idx(Ny+1,z)] = C[pos2idx(Ny,z)]

    #BC3: at y = Ny, convection only [accumulation of biomass at end]

    #y = Ny,  z = 0 .. Nz

  

    for j=0:Nz

        dM[pos2idx(Ny,j)] = ( + Dz * (M[pos2idx(Ny,max(0,j-1))] + M[pos2idx(Ny,min(Nz,j+1))] - 2*M[pos2idx(Ny,j)]) / dz(Ny,Temperature[pos2idx(Ny,j)],0)^2
                              - Vol(dz(Ny,Temperature[pos2idx(Ny,j)],0))*V_profile(dz(Ny,Temperature[pos2idx(Ny,j)],0)*j, Hght(Ny,Temperature[pos2idx(Ny,j)],0)) * (M[pos2idx(Ny,j)]/Vol(dz(Ny,Temperature[pos2idx(Ny,j)],0)) - M[pos2idx(Ny-1,j)]/Vol(dz(Ny-1,Temperature[pos2idx(Ny-1,j)],0))) / dy
                            )
    end

    #y = dy .. Ny-dy,  z = 0 or z = Nz
    for i=1:Ny-1
        

        dM[pos2idx(i,0)] =     ( Vol(dz(i,Temperature[pos2idx(i,0)],0))*Dy * (M[pos2idx(i-1,0)]/Vol(dz(i-1,Temperature[pos2idx(i-1,0)],0)) + M[pos2idx(i+1,0)]/Vol(dz(i+1,Temperature[pos2idx(i+1,0)],0)) - 2*M[pos2idx(i,0)]/Vol(dz(i,Temperature[pos2idx(i,0)],0))) / dy^2  #small
                               + Dz * (M[pos2idx(i,0)] + M[pos2idx(i,0+1)] - 2*M[pos2idx(i,0)]) /dz(i,Temperature[pos2idx(i,0)],0)^2
                               - Vol(dz(i,Temperature[pos2idx(i,0)],0))*V_profile(dz(i,Temperature[pos2idx(i,0)],0)*0, Hght(i,Temperature[pos2idx(i,0)],0)) * ( M[pos2idx(i,0)]/Vol(dz(i,Temperature[pos2idx(i,0)],0)) - M[pos2idx(i-1,0)]/Vol(dz(i-1,Temperature[pos2idx(i-1,0)],0)) )/dy
                               + mu(M[pos2idx(i,0:Nz)],dz(i,Temperature[pos2idx(i,0)],0), 0) * M[pos2idx(i,0)]
                               )

        dM[pos2idx(i,Nz)] =    ( Vol(dz(i,Temperature[pos2idx(i,Nz)],0))*Dy * (M[pos2idx(i-1,Nz)]/Vol(dz(i-1,Temperature[pos2idx(i-1,Nz)],0)) + M[pos2idx(i+1,Nz)]/Vol(dz(i+1,Temperature[pos2idx(i+1,Nz)],0)) - 2*M[pos2idx(i,Nz)]/Vol(dz(i,Temperature[pos2idx(i,Nz)],0))) / dy^2 #small
                               + Dz * (M[pos2idx(i,Nz-1)] + M[pos2idx(i,Nz)] - 2*M[pos2idx(i,Nz)]) / dz(i,Temperature[pos2idx(i,Nz)],0)^2
                               - Vol(dz(i,Temperature[pos2idx(i,Nz)],0))*V_profile(dz(i,Temperature[pos2idx(i,Nz)],0)*Nz, Hght(i,Temperature[pos2idx(i,Nz)],0)) * (M[pos2idx(i,Nz)]/Vol(dz(i,Temperature[pos2idx(i,Nz)],0)) - M[pos2idx(i-1,Nz)]/Vol(dz(i-1,Temperature[pos2idx(i-1,Nz)],0))) / dy
                               + mu(M[pos2idx(i,0:Nz)],dz(i,Temperature[pos2idx(i,Nz)],0), Nz) * M[pos2idx(i,Nz)]
                               )
    end

    for i=1:Ny-1
        for j=1:Nz-1

            dM[pos2idx(i,j)] = ( Vol(dz(i,Temperature[pos2idx(i,j)],0))*Dy * (M[pos2idx(i-1,j)]/Vol(dz(i-1,Temperature[pos2idx(i-1,j)],0)) + M[pos2idx(i+1,j)]/Vol(dz(i+1,Temperature[pos2idx(i+1,j)],0)) - 2*M[pos2idx(i,j)]/Vol(dz(i,Temperature[pos2idx(i,j)],0))) / dy^2
                               + Dz * (M[pos2idx(i,j-1)] + M[pos2idx(i,j+1)] - 2*M[pos2idx(i,j)]) / dz(i,Temperature[pos2idx(i,j)],0)^2
                               - Vol(dz(i,Temperature[pos2idx(i,j)],0))*V_profile(dz(i,Temperature[pos2idx(i,j)],0)*j, Hght(i,Temperature[pos2idx(i,j)],0)) * (M[pos2idx(i,j)]/Vol(dz(i,Temperature[pos2idx(i,j)],0)) - M[pos2idx(i-1,j)]/Vol(dz(i-1,Temperature[pos2idx(i-1,j)],0))) / dy
                               + mu(M[pos2idx(i,0:Nz)],dz(i,Temperature[pos2idx(i,j)],0), j) * M[pos2idx(i,j)]
                               )
        end
    end

    @views dX[1:Nelements] .= dM        # order matters! The @views operator takes a slice out of an array without making a copy.

    nothing
end
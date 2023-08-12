cd("C:/Users/slant/OneDrive/Desktop/Julia_2")
function PDE_AlgaeBiomass!(dX, M, Height, params, t)

    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data

    t_hour1 = floor(Int64, t)
    t_hour2 = floor(Int64, t)+1

    
    data_begin = params.data_begin

    GHI = GHI_Data[data_begin + t_hour1] * (t-t_hour1) + GHI_Data[data_begin + t_hour2] * (t_hour2 - t)
    Tamb = max.(((Tamb_Data[data_begin + t_hour1] * (t-t_hour1) + Tamb_Data[data_begin + t_hour2] * (t_hour2 - t))),0)


    Isat = params.max_biomass_light_saturation
    Cmax = params.max_biomass_concentration
    mu_max = params.max_biomass_specific_growth_rate
    Q = params.input_volumetric_flow_rate
    Vmax = params.input_max_flow_velocity
    V_prof = params.velocity_profile
    Pt_eng = params.change_potential_energy
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
    dz = zeros(Nelements,1)
    
    for i = 0:Ny
        for j = 0:Nz
            dz[pos2idx(i,j)] = Height[pos2idx(i,j)]/Nz
        end
    end



    ##convert to mass balance
    Vol(dz_y) = dy*dz_y*W
    Iave(Mz, dz_y, z) = GHI .* 0.45 .* (1 - min(sum(Mz[1:z]/(Vol(dz_y))).*dz_y./Cmax,1.0) ) #
    phiL(Mz, dz_y, z) = Iave(Mz, dz_y, z) .* exp(1 - Iave(Mz, dz_y, z)./Isat) ./ Isat
    mu(Mz, dz_y, z) = mu_max .* phiL(Mz, dz_y, z)
    V_profile(x,z,H) = V_prof(x,z,H)

    
        

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

        dM[pos2idx(Ny,j)] = ( + Dz * (M[pos2idx(Ny,max(0,j-1))] + M[pos2idx(Ny,min(Nz,j+1))] - 2*M[pos2idx(Ny,j)]) / dz[pos2idx(Ny,j)]^2
                              - Vol(dz[pos2idx(Ny,j)])*V_profile(Ny, min(j,Nz-1)*dz[pos2idx(Ny,j)], Height[pos2idx(Ny,j)]) * (M[pos2idx(Ny,j)]/Vol(dz[pos2idx(Ny,j)]) - M[pos2idx(Ny-1,j)]/Vol(dz[pos2idx(Ny-1,j)])) / dy
                            )
    end

    #y = dy .. Ny-dy,  z = 0 or z = Nz
    for i=1:Ny-1
        

        dM[pos2idx(i,0)] =     ( Vol(dz[pos2idx(i,0)])*Dy * (M[pos2idx(i-1,0)]/Vol(dz[pos2idx(i-1,0)]) + M[pos2idx(i+1,0)]/Vol(dz[pos2idx(i+1,0)]) - 2*M[pos2idx(i,0)]/Vol(dz[pos2idx(i,0)])) / dy^2  #small
                               + Dz * (M[pos2idx(i,0)] + M[pos2idx(i,0+1)] - 2*M[pos2idx(i,0)]) /dz[pos2idx(i,0)]^2
                               - Vol(dz[pos2idx(i,0)])*V_profile(i, 0, Height[pos2idx(i,0)]) * ( M[pos2idx(i,0)]/Vol(dz[pos2idx(i,0)]) - M[pos2idx(i-1,0)]/Vol(dz[pos2idx(i-1,0)]) )/dy
                               + mu(M[pos2idx(i,0:Nz)],dz[pos2idx(i,0)], 0) * M[pos2idx(i,0)]
                               )

        dM[pos2idx(i,Nz)] =    ( Vol(dz[pos2idx(i,Nz)])*Dy * (M[pos2idx(i-1,Nz)]/Vol(dz[pos2idx(i-1,Nz)]) + M[pos2idx(i+1,Nz)]/Vol(dz[pos2idx(i+1,Nz)]) - 2*M[pos2idx(i,Nz)]/Vol(dz[pos2idx(i,Nz)])) / dy^2 #small
                               + Dz * (M[pos2idx(i,Nz-1)] + M[pos2idx(i,Nz)] - 2*M[pos2idx(i,Nz)]) / dz[pos2idx(i,Nz)]^2
                               - Vol(dz[pos2idx(i,Nz)])*V_profile(i,Nz*dz[pos2idx(i,Nz)], Height[pos2idx(i,Nz)]) * (M[pos2idx(i,Nz)]/Vol(dz[pos2idx(i,Nz)]) - M[pos2idx(i-1,Nz)]/Vol(dz[pos2idx(i-1,Nz)])) / dy
                               + mu(M[pos2idx(i,0:Nz)],dz[pos2idx(i,Nz)], Nz) * M[pos2idx(i,Nz)]
                               )
    end

    for i=1:Ny-1
        for j=1:Nz-1

            dM[pos2idx(i,j)] = ( Vol(dz[pos2idx(i,j)])*Dy * (M[pos2idx(i-1,j)]/Vol(dz[pos2idx(i-1,j)]) + M[pos2idx(i+1,j)]/Vol(dz[pos2idx(i+1,j)]) - 2*M[pos2idx(i,j)]/Vol(dz[pos2idx(i,j)])) / dy^2
                               + Dz * (M[pos2idx(i,j-1)] + M[pos2idx(i,j+1)] - 2*M[pos2idx(i,j)]) / dz[pos2idx(i,j)]^2
                               - Vol(dz[pos2idx(i,j)])*V_profile(i,j*dz[pos2idx(i,j)], Height[pos2idx(i,j)]) * (M[pos2idx(i,j)]/Vol(dz[pos2idx(i,j)]) - M[pos2idx(i-1,j)]/Vol(dz[pos2idx(i-1,j)])) / dy
                               + mu(M[pos2idx(i,0:Nz)],dz[pos2idx(i,j)], j) * M[pos2idx(i,j)]
                               )
        end
    end

    @views dX[1:Nelements] .= dM        # order matters! The @views operator takes a slice out of an array without making a copy.

    nothing
end
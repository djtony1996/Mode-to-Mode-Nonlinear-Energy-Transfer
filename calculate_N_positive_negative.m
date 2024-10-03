function [N_positive, N_negative, N_positive_z, N_negative_z] = calculate_N_positive_negative(kx_0posi,ky_0posi,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,nz,dkx,dky,WEIGHT)

N_positive     = zeros(length(ky_0posi),length(kx_0posi));
N_negative     = zeros(length(ky_0posi),length(kx_0posi));
N_positive_z   = zeros(nz-1,length(ky_0posi),length(kx_0posi));   
N_negative_z   = zeros(nz-1,length(ky_0posi),length(kx_0posi));   

% calculate each nonlinear energy transfer
parfor kx_wave = 1: length(kx_0posi)
    disp(kx_wave)

    N_positive_temp     = zeros(length(ky_0posi),1);
    N_negative_temp     = zeros(length(ky_0posi),1);
    N_positive_z_temp   = zeros(nz-1,length(ky_0posi));   
    N_negative_z_temp   = zeros(nz-1,length(ky_0posi));   

    for ky_wave = 1: length(ky_0posi)
        k_x = kx_0posi(kx_wave);
        k_y = ky_0posi(ky_wave);

        if k_x == 0 && k_y == 0
           continue
        end

        if k_x == 0
            M_single = (-2) .* real(get_N_transfer_0b(k_y,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky,WEIGHT));
            M_single3= (-2) .* real(get_N_transfer_0b(k_y,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky));
        elseif k_y == 0
            M_single = (-2) .* real(get_N_transfer_a0(k_x,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky,WEIGHT));
            M_single3= (-2) .* real(get_N_transfer_a0(k_x,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky));
        else 
            M_single = (-2) .* real(get_N_transfer_ab(k_x,k_y,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky,WEIGHT));
            M_single3= (-2) .* real(get_N_transfer_ab(k_x,k_y,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky));
        end    
        
        N_positive_temp(ky_wave)     = sum(sum(M_single .* (M_single>0)));
        N_negative_temp(ky_wave)     = sum(sum(M_single .* (M_single<0)));
        N_positive_z_temp(:,ky_wave) = sum((M_single3.*(M_single3>0)),[2 3]);
        N_negative_z_temp(:,ky_wave) = sum((M_single3.*(M_single3<0)),[2 3]);
    end

    N_positive(:,kx_wave)     = N_positive_temp;
    N_negative(:,kx_wave)     = N_negative_temp;
    N_positive_z(:,:,kx_wave) = N_positive_z_temp;
    N_negative_z(:,:,kx_wave) = N_negative_z_temp;
end

N_positive = single(N_positive);
N_negative = single(N_negative);
N_positive_z = single(N_positive_z);
N_negative_z = single(N_negative_z);

end

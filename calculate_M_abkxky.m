function [M_4d] = calculate_M_abkxky(kx_0posi,ky_0posi,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky,WEIGHT)

M_4d     = zeros(length(ky_0posi),length(kx_0posi),length(ky_0posi),length(kx_0posi));
% M_5d     = zeros(nz-1,length(ky_0posi),length(kx_0posi),length(ky_0posi),length(kx_0posi));   

% calculate each nonlinear energy transfer
parfor kx_wave = 1: length(kx_0posi)
    disp(kx_wave)
    M_3d_temp = zeros(length(ky_0posi),length(kx_0posi),length(ky_0posi));
    for ky_wave = 1: length(ky_0posi)
        k_x = kx_0posi(kx_wave);
        k_y = ky_0posi(ky_wave);

        if k_x == 0 && k_y == 0
           continue
        end

        if k_x == 0
            M_single = (-2) .* real(get_N_transfer_0b(k_y,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky,WEIGHT));
            % M_single3= (-2) .* real(get_N_transfer_0b(k_y,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky));
        elseif k_y == 0
            M_single = (-2) .* real(get_N_transfer_a0(k_x,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky,WEIGHT));
            % M_single3= (-2) .* real(get_N_transfer_a0(k_x,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky));
        else 
            M_single = (-2) .* real(get_N_transfer_ab(k_x,k_y,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky,WEIGHT));
            % M_single3= (-2) .* real(get_N_transfer_ab(k_x,k_y,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky));
        end    
        % M_matrix(:,ky_wave-1+(kx_wave-1)*length(ky_0posi)) = M_single(2:end);  
        M_3d_temp(:,:,ky_wave)    = M_single;
        % M_5d(:,:,:,ky_wave,kx_wave) = M_single3;
    end
    M_4d(:,:,:,kx_wave) = M_3d_temp;
end

end
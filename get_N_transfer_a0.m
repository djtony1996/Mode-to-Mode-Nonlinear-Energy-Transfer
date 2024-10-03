function [N_transfer] = get_N_transfer_a0(k_x,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,Nx,Ny,dkx,dky,WEIGHT)

N_transfer   = zeros(size(u_F,1),length(ky_posi)+1,length(kx_posi)+1);

[index_3_x, index_3_y]     = get_index(-k_x,0,Nx,Ny,dkx,dky);

for index_k_a = 1: length(kx_posi)
    for index_k_b = 1: length(ky_posi)
        k_a = kx_posi(index_k_a);
        k_b = ky_posi(index_k_b);

        [index_1_1_x, index_1_1_y] = get_index(k_x-k_a,-k_b,Nx,Ny,dkx,dky);
        [index_1_2_x, index_1_2_y] = get_index(k_x+k_a,+k_b,Nx,Ny,dkx,dky);
        [index_1_3_x, index_1_3_y] = get_index(k_x+k_a,-k_b,Nx,Ny,dkx,dky);
        [index_1_4_x, index_1_4_y] = get_index(k_x-k_a,+k_b,Nx,Ny,dkx,dky);        

        [index_2_1_x, index_2_1_y] = get_index(+k_a,+k_b ,Nx,Ny,dkx,dky);
        [index_2_2_x, index_2_2_y] = get_index(-k_a,-k_b,Nx,Ny,dkx,dky);
        [index_2_3_x, index_2_3_y] = get_index(-k_a,+k_b ,Nx,Ny,dkx,dky);
        [index_2_4_x, index_2_4_y] = get_index(+k_a,-k_b,Nx,Ny,dkx,dky);
   
        [N_transfer(:,round(k_b/dky+1),round(k_a/dkx+1))] = get_uj_duidxj_ui(u_F(:,index_1_1_y,index_1_1_x),v_F(:,index_1_1_y,index_1_1_x),w_F(:,index_1_1_y,index_1_1_x),duF_dx(:,index_2_1_y,index_2_1_x),duF_dy(:,index_2_1_y,index_2_1_x),duF_dz(:,index_2_1_y,index_2_1_x),dvF_dx(:,index_2_1_y,index_2_1_x),dvF_dy(:,index_2_1_y,index_2_1_x),dvF_dz(:,index_2_1_y,index_2_1_x),dwF_dx(:,index_2_1_y,index_2_1_x),dwF_dy(:,index_2_1_y,index_2_1_x),dwF_dz(:,index_2_1_y,index_2_1_x),u_F(:,index_3_y,index_3_x),v_F(:,index_3_y,index_3_x),w_F(:,index_3_y,index_3_x))...
                                                          + get_uj_duidxj_ui(u_F(:,index_1_2_y,index_1_2_x),v_F(:,index_1_2_y,index_1_2_x),w_F(:,index_1_2_y,index_1_2_x),duF_dx(:,index_2_2_y,index_2_2_x),duF_dy(:,index_2_2_y,index_2_2_x),duF_dz(:,index_2_2_y,index_2_2_x),dvF_dx(:,index_2_2_y,index_2_2_x),dvF_dy(:,index_2_2_y,index_2_2_x),dvF_dz(:,index_2_2_y,index_2_2_x),dwF_dx(:,index_2_2_y,index_2_2_x),dwF_dy(:,index_2_2_y,index_2_2_x),dwF_dz(:,index_2_2_y,index_2_2_x),u_F(:,index_3_y,index_3_x),v_F(:,index_3_y,index_3_x),w_F(:,index_3_y,index_3_x))...
                                                          + get_uj_duidxj_ui(u_F(:,index_1_3_y,index_1_3_x),v_F(:,index_1_3_y,index_1_3_x),w_F(:,index_1_3_y,index_1_3_x),duF_dx(:,index_2_3_y,index_2_3_x),duF_dy(:,index_2_3_y,index_2_3_x),duF_dz(:,index_2_3_y,index_2_3_x),dvF_dx(:,index_2_3_y,index_2_3_x),dvF_dy(:,index_2_3_y,index_2_3_x),dvF_dz(:,index_2_3_y,index_2_3_x),dwF_dx(:,index_2_3_y,index_2_3_x),dwF_dy(:,index_2_3_y,index_2_3_x),dwF_dz(:,index_2_3_y,index_2_3_x),u_F(:,index_3_y,index_3_x),v_F(:,index_3_y,index_3_x),w_F(:,index_3_y,index_3_x))...
                                                          + get_uj_duidxj_ui(u_F(:,index_1_4_y,index_1_4_x),v_F(:,index_1_4_y,index_1_4_x),w_F(:,index_1_4_y,index_1_4_x),duF_dx(:,index_2_4_y,index_2_4_x),duF_dy(:,index_2_4_y,index_2_4_x),duF_dz(:,index_2_4_y,index_2_4_x),dvF_dx(:,index_2_4_y,index_2_4_x),dvF_dy(:,index_2_4_y,index_2_4_x),dvF_dz(:,index_2_4_y,index_2_4_x),dwF_dx(:,index_2_4_y,index_2_4_x),dwF_dy(:,index_2_4_y,index_2_4_x),dwF_dz(:,index_2_4_y,index_2_4_x),u_F(:,index_3_y,index_3_x),v_F(:,index_3_y,index_3_x),w_F(:,index_3_y,index_3_x));

%         [temp1,temp1_x,temp1_y,temp1_z] = get_uj_duidxj_ui(u_F(:,index_1_1_y,index_1_1_x),v_F(:,index_1_1_y,index_1_1_x),w_F(:,index_1_1_y,index_1_1_x),duF_dx(:,index_2_1_y,index_2_1_x),duF_dy(:,index_2_1_y,index_2_1_x),duF_dz(:,index_2_1_y,index_2_1_x),dvF_dx(:,index_2_1_y,index_2_1_x),dvF_dy(:,index_2_1_y,index_2_1_x),dvF_dz(:,index_2_1_y,index_2_1_x),dwF_dx(:,index_2_1_y,index_2_1_x),dwF_dy(:,index_2_1_y,index_2_1_x),dwF_dz(:,index_2_1_y,index_2_1_x),u_F(:,index_3_y,index_3_x),v_F(:,index_3_y,index_3_x),w_F(:,index_3_y,index_3_x));
%         [temp2,temp2_x,temp2_y,temp2_z] = get_uj_duidxj_ui(u_F(:,index_1_2_y,index_1_2_x),v_F(:,index_1_2_y,index_1_2_x),w_F(:,index_1_2_y,index_1_2_x),duF_dx(:,index_2_2_y,index_2_2_x),duF_dy(:,index_2_2_y,index_2_2_x),duF_dz(:,index_2_2_y,index_2_2_x),dvF_dx(:,index_2_2_y,index_2_2_x),dvF_dy(:,index_2_2_y,index_2_2_x),dvF_dz(:,index_2_2_y,index_2_2_x),dwF_dx(:,index_2_2_y,index_2_2_x),dwF_dy(:,index_2_2_y,index_2_2_x),dwF_dz(:,index_2_2_y,index_2_2_x),u_F(:,index_3_y,index_3_x),v_F(:,index_3_y,index_3_x),w_F(:,index_3_y,index_3_x));
%         [temp3,temp3_x,temp3_y,temp3_z] = get_uj_duidxj_ui(u_F(:,index_1_3_y,index_1_3_x),v_F(:,index_1_3_y,index_1_3_x),w_F(:,index_1_3_y,index_1_3_x),duF_dx(:,index_2_3_y,index_2_3_x),duF_dy(:,index_2_3_y,index_2_3_x),duF_dz(:,index_2_3_y,index_2_3_x),dvF_dx(:,index_2_3_y,index_2_3_x),dvF_dy(:,index_2_3_y,index_2_3_x),dvF_dz(:,index_2_3_y,index_2_3_x),dwF_dx(:,index_2_3_y,index_2_3_x),dwF_dy(:,index_2_3_y,index_2_3_x),dwF_dz(:,index_2_3_y,index_2_3_x),u_F(:,index_3_y,index_3_x),v_F(:,index_3_y,index_3_x),w_F(:,index_3_y,index_3_x));
%         [temp4,temp4_x,temp4_y,temp4_z] = get_uj_duidxj_ui(u_F(:,index_1_4_y,index_1_4_x),v_F(:,index_1_4_y,index_1_4_x),w_F(:,index_1_4_y,index_1_4_x),duF_dx(:,index_2_4_y,index_2_4_x),duF_dy(:,index_2_4_y,index_2_4_x),duF_dz(:,index_2_4_y,index_2_4_x),dvF_dx(:,index_2_4_y,index_2_4_x),dvF_dy(:,index_2_4_y,index_2_4_x),dvF_dz(:,index_2_4_y,index_2_4_x),dwF_dx(:,index_2_4_y,index_2_4_x),dwF_dy(:,index_2_4_y,index_2_4_x),dwF_dz(:,index_2_4_y,index_2_4_x),u_F(:,index_3_y,index_3_x),v_F(:,index_3_y,index_3_x),w_F(:,index_3_y,index_3_x));
%         N_transfer(:,k_b/dky+1,k_a/dkx+1)  = temp1 + temp2 + temp3 + temp4;
%         N_transfer_x(:,k_b/dky+1,k_a/dkx+1)= temp1_x + temp2_x + temp3_x + temp4_x;
%         N_transfer_y(:,k_b/dky+1,k_a/dkx+1)= temp1_y + temp2_y + temp3_y + temp4_y;
%         N_transfer_z(:,k_b/dky+1,k_a/dkx+1)= temp1_z + temp2_z + temp3_z + temp4_z;
    end
end

for index_k_a = 1: length(kx_posi)
    k_a = kx_posi(index_k_a);

    [index_1_1_x, index_1_1_y] = get_index(k_x-k_a,0,Nx,Ny,dkx,dky);
    [index_1_2_x, index_1_2_y] = get_index(k_x+k_a,0,Nx,Ny,dkx,dky);     

    [index_2_1_x, index_2_1_y] = get_index(+k_a,0,Nx,Ny,dkx,dky);
    [index_2_2_x, index_2_2_y] = get_index(-k_a,0,Nx,Ny,dkx,dky);

    [N_transfer(:,1,round(k_a/dkx+1))] = get_uj_duidxj_ui(u_F(:,index_1_1_y,index_1_1_x),v_F(:,index_1_1_y,index_1_1_x),w_F(:,index_1_1_y,index_1_1_x),duF_dx(:,index_2_1_y,index_2_1_x),duF_dy(:,index_2_1_y,index_2_1_x),duF_dz(:,index_2_1_y,index_2_1_x),dvF_dx(:,index_2_1_y,index_2_1_x),dvF_dy(:,index_2_1_y,index_2_1_x),dvF_dz(:,index_2_1_y,index_2_1_x),dwF_dx(:,index_2_1_y,index_2_1_x),dwF_dy(:,index_2_1_y,index_2_1_x),dwF_dz(:,index_2_1_y,index_2_1_x),u_F(:,index_3_y,index_3_x),v_F(:,index_3_y,index_3_x),w_F(:,index_3_y,index_3_x))...
                                       + get_uj_duidxj_ui(u_F(:,index_1_2_y,index_1_2_x),v_F(:,index_1_2_y,index_1_2_x),w_F(:,index_1_2_y,index_1_2_x),duF_dx(:,index_2_2_y,index_2_2_x),duF_dy(:,index_2_2_y,index_2_2_x),duF_dz(:,index_2_2_y,index_2_2_x),dvF_dx(:,index_2_2_y,index_2_2_x),dvF_dy(:,index_2_2_y,index_2_2_x),dvF_dz(:,index_2_2_y,index_2_2_x),dwF_dx(:,index_2_2_y,index_2_2_x),dwF_dy(:,index_2_2_y,index_2_2_x),dwF_dz(:,index_2_2_y,index_2_2_x),u_F(:,index_3_y,index_3_x),v_F(:,index_3_y,index_3_x),w_F(:,index_3_y,index_3_x));

end

for index_k_b = 1: length(ky_posi)
    k_b = ky_posi(index_k_b);

    [index_1_1_x, index_1_1_y] = get_index(k_x,-k_b,Nx,Ny,dkx,dky);
    [index_1_2_x, index_1_2_y] = get_index(k_x,+k_b,Nx,Ny,dkx,dky);

    [index_2_1_x, index_2_1_y] = get_index(0,k_b,Nx,Ny,dkx,dky);
    [index_2_2_x, index_2_2_y] = get_index(0,-k_b,Nx,Ny,dkx,dky);

    [N_transfer(:,round(k_b/dky+1),1)] = get_uj_duidxj_ui(u_F(:,index_1_1_y,index_1_1_x),v_F(:,index_1_1_y,index_1_1_x),w_F(:,index_1_1_y,index_1_1_x),duF_dx(:,index_2_1_y,index_2_1_x),duF_dy(:,index_2_1_y,index_2_1_x),duF_dz(:,index_2_1_y,index_2_1_x),dvF_dx(:,index_2_1_y,index_2_1_x),dvF_dy(:,index_2_1_y,index_2_1_x),dvF_dz(:,index_2_1_y,index_2_1_x),dwF_dx(:,index_2_1_y,index_2_1_x),dwF_dy(:,index_2_1_y,index_2_1_x),dwF_dz(:,index_2_1_y,index_2_1_x),u_F(:,index_3_y,index_3_x),v_F(:,index_3_y,index_3_x),w_F(:,index_3_y,index_3_x))...
                                       + get_uj_duidxj_ui(u_F(:,index_1_2_y,index_1_2_x),v_F(:,index_1_2_y,index_1_2_x),w_F(:,index_1_2_y,index_1_2_x),duF_dx(:,index_2_2_y,index_2_2_x),duF_dy(:,index_2_2_y,index_2_2_x),duF_dz(:,index_2_2_y,index_2_2_x),dvF_dx(:,index_2_2_y,index_2_2_x),dvF_dy(:,index_2_2_y,index_2_2_x),dvF_dz(:,index_2_2_y,index_2_2_x),dwF_dx(:,index_2_2_y,index_2_2_x),dwF_dy(:,index_2_2_y,index_2_2_x),dwF_dz(:,index_2_2_y,index_2_2_x),u_F(:,index_3_y,index_3_x),v_F(:,index_3_y,index_3_x),w_F(:,index_3_y,index_3_x));

end

if nargin == 20
    N_transfer = squeeze(pagemtimes(WEIGHT,N_transfer));
end

end
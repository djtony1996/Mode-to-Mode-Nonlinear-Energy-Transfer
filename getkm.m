% prepare the matrices for calculating the x,y derivatives using spectral
% method
function [kx_m,ky_m] = getkm(alpha_array,beta_array,nx,ny,dkx,dky)
    kx_m = zeros(ny,nx);
    ky_m = zeros(ny,nx);
    
    for k_array = 1: length(alpha_array)

        index_kx = round(alpha_array(k_array)/dkx) + 1;
        index_ky = round(beta_array(k_array)/dky) + 1;

        if alpha_array(k_array)<0
            kx_m(index_ky,nx+index_kx) = alpha_array(k_array);
            ky_m(index_ky,nx+index_kx) = beta_array(k_array);
        else
            kx_m(index_ky,index_kx) = alpha_array(k_array);
            ky_m(index_ky,index_kx) = beta_array(k_array);
        end
    end
    kx_m(ny:-1:ny/2+2,:) = kx_m(2:1:ny/2,:);
    ky_m(ny:-1:ny/2+2,:) = -ky_m(2:1:ny/2,:);
end
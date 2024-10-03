% note that ky_array only contains positive wavenumbers
function [F_vector] = getF(F_matrix, kx_array, ky_array, dkx, dky)
    F_vector = zeros(length(kx_array),1);
    
    for k_array = 1: length(kx_array)

        index_kx = round(kx_array(k_array)/dkx) + 1;
        index_ky = round(ky_array(k_array)/dky) + 1;

        if kx_array(k_array)<0
            F_vector(k_array) = F_matrix(index_ky,size(F_matrix,2)+index_kx); 
        else
            F_vector(k_array) = F_matrix(index_ky,index_kx); 
        end
    end
end
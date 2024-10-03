function [M_4d_avg_restrict] = get_restrict_source_recipient_M(M_4d_avg,kx_source_restriction_array,ky_source_restriction_array,kx_recipient_restriction_array,ky_recipient_restriction_array)

M_4d_avg_restrict = M_4d_avg;

% restrict source
for k_index = 1: length(kx_source_restriction_array) 
    M_single = M_4d_avg_restrict(:,:,ky_source_restriction_array(k_index),kx_source_restriction_array(k_index));
    M_single(M_single<0) = 0;
    M_4d_avg_restrict(:,:,ky_source_restriction_array(k_index),kx_source_restriction_array(k_index)) = M_single;
end

% delete M representing receiving the discarded nonlinear energy
filter = zeros(size(M_4d_avg,1),size(M_4d_avg,2));
ind = sub2ind(size(filter), kx_source_restriction_array, ky_source_restriction_array);
filter(ind) = 1;
for kx_index = 1: size(M_4d_avg,4)
    for ky_index = 1: size(M_4d_avg,3)
        M_single = M_4d_avg_restrict(:,:,ky_index,kx_index);
        filter1 = filter .* (M_single>0);
        M_single = M_single .* (~filter1);
        M_4d_avg_restrict(:,:,ky_index,kx_index) = M_single;
    end
end

% restrict recipient
for k_index = 1: length(kx_recipient_restriction_array) 
    M_single = M_4d_avg_restrict(:,:,ky_recipient_restriction_array(k_index),kx_recipient_restriction_array(k_index));
    M_single(M_single>0) = 0;
    M_4d_avg_restrict(:,:,ky_recipient_restriction_array(k_index),kx_recipient_restriction_array(k_index)) = M_single;
end

% delete M representing losing the discarded nonlinear energy
filter = zeros(size(M_4d_avg,1),size(M_4d_avg,2));
ind = sub2ind(size(filter), kx_recipient_restriction_array, ky_recipient_restriction_array);
filter(ind) = 1;
for kx_index = 1: size(M_4d_avg,4)
    for ky_index = 1: size(M_4d_avg,3)
        M_single = M_4d_avg_restrict(:,:,ky_index,kx_index);
        filter1 = filter .* (M_single<0);
        M_single = M_single .* (~filter1);
        M_4d_avg_restrict(:,:,ky_index,kx_index) = M_single;
    end
end

end
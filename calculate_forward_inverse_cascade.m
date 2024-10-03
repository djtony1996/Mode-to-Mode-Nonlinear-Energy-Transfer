function calculate_forward_inverse_cascade(Retau,nc_array,kx_max,ky_max)

load(['full',num2str(Retau),'_mean.mat'],'dkx','dky'); 
load(['M_',num2str(Retau),'_avg_final.mat'],'M_4d_avg_final'); M_4d_avg = M_4d_avg_final; clear M_4d_avg_final

N_forward_cascade_sum_positive    = zeros(length(nc_array),1); 
N_forward_cascade_sum_negative    = zeros(length(nc_array),1); 
N_inverse_cascade_sum_positive    = zeros(length(nc_array),1); 
N_inverse_cascade_sum_negative    = zeros(length(nc_array),1); 

for k_index = 1: length(nc_array)
k_index
nc = nc_array(k_index);
kxc = nc;
kyc = 2*nc;

% for forward cascade
[kx_recipient_restriction_array,ky_recipient_restriction_array] = get_rectangular_restriction(kxc,kyc,dkx,dky);
[kx_source_restriction_array,ky_source_restriction_array] = get_7_restriction(kxc,kx_max,kyc,ky_max,dkx,dky);

[M_4d_avg_restrict] = get_restrict_source_recipient_M(M_4d_avg,kx_source_restriction_array,ky_source_restriction_array,kx_recipient_restriction_array,ky_recipient_restriction_array);

N_positive_restrict = zeros(size(M_4d_avg_restrict,1),size(M_4d_avg_restrict,2));
N_negative_restrict = zeros(size(M_4d_avg_restrict,1),size(M_4d_avg_restrict,2));

for kx_index = 1: size(M_4d_avg,4)
    for ky_index = 1: size(M_4d_avg,3)
        M_single2 = M_4d_avg_restrict(:,:,ky_index,kx_index);
        N_positive_restrict(ky_index,kx_index) = sum(sum(M_single2 .* (M_single2>0)));
        N_negative_restrict(ky_index,kx_index) = sum(sum(M_single2 .* (M_single2<0)));
    end
end

N_forward_cascade_sum_positive(k_index) = sum(sum(N_positive_restrict));
N_forward_cascade_sum_negative(k_index) = sum(sum(N_negative_restrict));




% for inverse cascade
[kx_source_restriction_array,ky_source_restriction_array] = get_rectangular_restriction(kxc,kyc,dkx,dky);
[kx_recipient_restriction_array,ky_recipient_restriction_array] = get_7_restriction(kxc,kx_max,kyc,ky_max,dkx,dky);

[M_4d_avg_restrict] = get_restrict_source_recipient_M(M_4d_avg,kx_source_restriction_array,ky_source_restriction_array,kx_recipient_restriction_array,ky_recipient_restriction_array);

N_positive_restrict = zeros(size(M_4d_avg_restrict,1),size(M_4d_avg_restrict,2));
N_negative_restrict = zeros(size(M_4d_avg_restrict,1),size(M_4d_avg_restrict,2));

for kx_index = 1: size(M_4d_avg,4)
    for ky_index = 1: size(M_4d_avg,3)
        M_single2 = M_4d_avg_restrict(:,:,ky_index,kx_index);
        N_positive_restrict(ky_index,kx_index) = sum(sum(M_single2 .* (M_single2>0)));
        N_negative_restrict(ky_index,kx_index) = sum(sum(M_single2 .* (M_single2<0)));
    end
end

N_inverse_cascade_sum_positive(k_index) = sum(sum(N_positive_restrict));
N_inverse_cascade_sum_negative(k_index) = sum(sum(N_negative_restrict));
end

N_positive = zeros(size(M_4d_avg_restrict,1),size(M_4d_avg_restrict,2));
N_negative = zeros(size(M_4d_avg_restrict,1),size(M_4d_avg_restrict,2));
for kx_index = 1: size(M_4d_avg,4)
    for ky_index = 1: size(M_4d_avg,3)
        M_single1 = M_4d_avg(:,:,ky_index,kx_index);
        N_positive(ky_index,kx_index) = sum(sum(M_single1 .* (M_single1>0)));
        N_negative(ky_index,kx_index) = sum(sum(M_single1 .* (M_single1<0)));
    end
end

N_sum_positive = sum(sum(N_positive));
N_sum_negative = sum(sum(N_negative));


savename = ['N_forward_inverse_',num2str(Retau),'.mat'];
save(savename,'N_forward_cascade_sum_positive','N_forward_cascade_sum_negative','N_inverse_cascade_sum_positive','N_inverse_cascade_sum_negative','N_sum_positive','N_sum_negative','nc_array','-v7.3')





clear,clc

fft_x = @(x) fft(x, [], 3) / size(x, 3);
fft_y = @(x) fft(x, [], 2) / size(x, 2);
fft_xy = @(x) fft_x(fft_y(x));

ifft_x = @(x) ifft(x, [], 3) * size(x, 3);
ifft_y = @(x) ifft(x, [], 2) * size(x, 2);
ifft_xy = @(x) ifft_x(ifft_y(x));

Retau = 590;
load(strcat('full',num2str(Retau),'_mean.mat'));
u_mean = channelRe.Up(2:end-1);
dUdz = channelRe.Up_diff1(2:end-1);

if Retau == 180
    dirname = 'grid_180/112x112x150';
else
    dirname = 'grid_590/384x384x500';
end

[Diff,zc] = cheb(nz);
Diff = Diff(2:end-1,2:end-1);
% [~,WEIGHT] = clenCurt(nz);
% WEIGHT = WEIGHT(2:end-1);

k_edge = 190;
kx = [0:dkx:(k_edge*dkx),-(k_edge*dkx):dkx:-dkx];
ky = [0:dky:(k_edge*dky)];
kx_array = repelem(kx,length(ky));
ky_array = repmat(ky,1,length(kx));
kx_array = kx_array.';
ky_array = ky_array.';
[kx_m,ky_m] = getkm(kx_array,ky_array,nx,ny,dkx,dky);

% read_array = [100000:5000:230000];
read_array = [100000:50000:1000000];

k_wave_limit_array = [5:5:190];
Dissp_z_avg= zeros(nz-1,length(k_wave_limit_array));

for k_array = 1:length(read_array)
    k_array
    u = zeros(nzDNS+2,ny,nx);
    v = zeros(nzDNS+2,ny,nx);
    w = zeros(nzDNS+1,ny,nx);
    loadname = strcat('../ChanFast/',dirname,'/outputdir/u_it',num2str(read_array(k_array),'%.0f'),'.dat');
    u(2:end-1,:,:) = read_bin(loadname, [nzDNS, ny, nx]);
    u = interp1(zp,u,zc);
    % u = u - mean(mean(u,3),2);
    u = u(2:end-1,:,:);
    u = u - u_mean;
    u = permute(u,[3 2 1]);
    u = interp1(xu,u,xp,'linear','extrap');
    u = permute(u,[3 2 1]);
    loadname = strcat('../ChanFast/',dirname,'/outputdir/v_it',num2str(read_array(k_array),'%.0f'),'.dat');
    v(2:end-1,:,:) = read_bin(loadname, [nzDNS, ny, nx]);
    v = interp1(zp,v,zc);
    % v = v - mean(mean(v,3),2);
    v = v(2:end-1,:,:);
    v = permute(v,[2 1 3]);
    v = interp1(yv,v,yp,'linear','extrap');
    v = permute(v, [2 1 3]);
    loadname = strcat('../ChanFast/',dirname,'/outputdir/w_it',num2str(read_array(k_array),'%.0f'),'.dat');
    w(2:end,:,:) = read_bin(loadname, [nzDNS, ny, nx]);
    w = interp1(zw,w,zc);
    % w = w - mean(mean(w,3),2);
    w = w(2:end-1,:,:);

    u_F  = fft_xy(u);
    v_F  = fft_xy(v);
    w_F  = fft_xy(w);
    
    Dissp_z = zeros(nz-1,length(k_wave_limit_array));
    for k_wave_limit_index = 1:length(k_wave_limit_array)
        k_wave_limit = k_wave_limit_array(k_wave_limit_index);

        [u_F_filtered] = get_filtered_velocity(u_F,k_wave_limit);
        [v_F_filtered] = get_filtered_velocity(v_F,k_wave_limit);
        [w_F_filtered] = get_filtered_velocity(w_F,k_wave_limit);
       
        [duF_dx_filtered,duF_dy_filtered,duF_dz_filtered]    = get_3d(u_F_filtered,Diff,kx_m,ky_m);
        [dvF_dx_filtered,dvF_dy_filtered,dvF_dz_filtered]    = get_3d(v_F_filtered,Diff,kx_m,ky_m);
        [dwF_dx_filtered,dwF_dy_filtered,dwF_dz_filtered]    = get_3d(w_F_filtered,Diff,kx_m,ky_m); 
       
        % calculate the energy in physical space
        Dissxp_3= -(ifft_xy(duF_dx_filtered).*ifft_xy(duF_dx_filtered)+ifft_xy(duF_dy_filtered).*ifft_xy(duF_dy_filtered)+ifft_xy(duF_dz_filtered).*ifft_xy(duF_dz_filtered))./Retau;
        Dissyp_3= -(ifft_xy(dvF_dx_filtered).*ifft_xy(dvF_dx_filtered)+ifft_xy(dvF_dy_filtered).*ifft_xy(dvF_dy_filtered)+ifft_xy(dvF_dz_filtered).*ifft_xy(dvF_dz_filtered))./Retau;
        Disszp_3= -(ifft_xy(dwF_dx_filtered).*ifft_xy(dwF_dx_filtered)+ifft_xy(dwF_dy_filtered).*ifft_xy(dwF_dy_filtered)+ifft_xy(dwF_dz_filtered).*ifft_xy(dwF_dz_filtered))./Retau;
        Dissp_3 = Dissxp_3 + Dissyp_3 + Disszp_3;
        Dissp_z(:,k_wave_limit_index)= mean(Dissp_3,[2 3]);
    end
    

    % calculate the average quantities
    Dissp_z_avg = (k_array-1).*Dissp_z_avg./k_array + Dissp_z./k_array;
end


savename = strcat('dissipation_LES',num2str(Retau),'.mat');
save(savename,'Dissp_z_avg','-v7.3')

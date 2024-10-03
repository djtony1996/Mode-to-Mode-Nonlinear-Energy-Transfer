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
[~,WEIGHT] = clenCurt(nz);
WEIGHT = WEIGHT(2:end-1);

k_edge = 150;
kx = [0:dkx:(k_edge*dkx),-(k_edge*dkx):dkx:-dkx];
ky = [0:dky:(k_edge*dky)];
kx_array = repelem(kx,length(ky));
ky_array = repmat(ky,1,length(kx));
kx_array = kx_array.';
ky_array = ky_array.';

read_array = [100000:5000:230000];

ProF_avg = zeros(k_edge,k_edge);
DissF_avg= zeros(k_edge,k_edge);
NonTF_avg= zeros(k_edge,k_edge);
ProF_z_avg = zeros(nz-1,k_edge,k_edge);
DissF_z_avg= zeros(nz-1,k_edge,k_edge);
NonTF_z_avg= zeros(nz-1,k_edge,k_edge);
Prop_z_avg = zeros(nz-1,1);
Dissp_z_avg= zeros(nz-1,1);
NonTp_z_avg= zeros(nz-1,1); 

[kx_m,ky_m] = getkm(kx_array,ky_array,nx,ny,dkx,dky);

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
   
    [duF_dx,duF_dy,duF_dz]    = get_3d(u_F,Diff,kx_m,ky_m);
    [dvF_dx,dvF_dy,dvF_dz]    = get_3d(v_F,Diff,kx_m,ky_m);
    [dwF_dx,dwF_dy,dwF_dz]    = get_3d(w_F,Diff,kx_m,ky_m); 
   
    % calculate the energy in physical space
    Prop_3  = -u.*w.*dUdz;
    Dissxp_3= -(ifft_xy(duF_dx).*ifft_xy(duF_dx)+ifft_xy(duF_dy).*ifft_xy(duF_dy)+ifft_xy(duF_dz).*ifft_xy(duF_dz))./Retau;
    Dissyp_3= -(ifft_xy(dvF_dx).*ifft_xy(dvF_dx)+ifft_xy(dvF_dy).*ifft_xy(dvF_dy)+ifft_xy(dvF_dz).*ifft_xy(dvF_dz))./Retau;
    Disszp_3= -(ifft_xy(dwF_dx).*ifft_xy(dwF_dx)+ifft_xy(dwF_dy).*ifft_xy(dwF_dy)+ifft_xy(dwF_dz).*ifft_xy(dwF_dz))./Retau;
    Dissp_3 = Dissxp_3 + Dissyp_3 + Disszp_3;
    NonTxp_3= -(u.*u.*ifft_xy(duF_dx) + u.*v.*ifft_xy(duF_dy) + u.*w.*ifft_xy(duF_dz));
    NonTyp_3= -(v.*u.*ifft_xy(dvF_dx) + v.*v.*ifft_xy(dvF_dy) + v.*w.*ifft_xy(dvF_dz));
    NonTzp_3= -(w.*u.*ifft_xy(dwF_dx) + w.*v.*ifft_xy(dwF_dy) + w.*w.*ifft_xy(dwF_dz));
    NonTp_3 = NonTxp_3 + NonTyp_3 + NonTzp_3;
    
    Prop_z = mean(Prop_3,[2 3]);
    Dissp_z= mean(Dissp_3,[2 3]);
    NonTp_z= mean(NonTp_3,[2 3]);

    % calculate the energy in Fourier space
    [ProF, DissF, NonTF, ProF_z, DissF_z, NonTF_z] = get_three_energy(u,v,w,kx_array,ky_array,nx,ny,nz,dkx,dky,Retau,dUdz,Diff);
    NonTF = real(NonTF); NonTF_z = real(NonTF_z);

    % calculate energy in positive kx and ky
    ProF(2:k_edge,2:k_edge) = ProF(2:k_edge,2:k_edge) + ProF(2:k_edge,(end:-1:end-k_edge+2));
    ProF = ProF(1:k_edge,1:k_edge);
    DissF(2:k_edge,2:k_edge) = DissF(2:k_edge,2:k_edge) + DissF(2:k_edge,(end:-1:end-k_edge+2));
    DissF = DissF(1:k_edge,1:k_edge);
    NonTF(2:k_edge,2:k_edge) = NonTF(2:k_edge,2:k_edge) + NonTF(2:k_edge,(end:-1:end-k_edge+2));
    NonTF = NonTF(1:k_edge,1:k_edge);
   
    ProF_z(:,2:k_edge,2:k_edge) = ProF_z(:,2:k_edge,2:k_edge) + ProF_z(:,2:k_edge,(end:-1:end-k_edge+2));
    ProF_z = ProF_z(:,1:k_edge,1:k_edge);
    DissF_z(:,2:k_edge,2:k_edge) = DissF_z(:,2:k_edge,2:k_edge) + DissF_z(:,2:k_edge,(end:-1:end-k_edge+2));
    DissF_z = DissF_z(:,1:k_edge,1:k_edge);
    NonTF_z(:,2:k_edge,2:k_edge) = NonTF_z(:,2:k_edge,2:k_edge) + NonTF_z(:,2:k_edge,(end:-1:end-k_edge+2));
    NonTF_z = NonTF_z(:,1:k_edge,1:k_edge);

    % calculate the average quantities
    Prop_z_avg  = (k_array-1).*Prop_z_avg./k_array + Prop_z./k_array;
    Dissp_z_avg = (k_array-1).*Dissp_z_avg./k_array + Dissp_z./k_array;
    NonTp_z_avg = (k_array-1).*NonTp_z_avg./k_array + NonTp_z./k_array;
    ProF_avg  = (k_array-1).*ProF_avg./k_array + ProF./k_array;
    DissF_avg = (k_array-1).*DissF_avg./k_array + DissF./k_array;
    NonTF_avg = (k_array-1).*NonTF_avg./k_array + NonTF./k_array;
    ProF_z_avg  = (k_array-1).*ProF_z_avg./k_array + ProF_z./k_array;
    DissF_z_avg = (k_array-1).*DissF_z_avg./k_array + DissF_z./k_array;
    NonTF_z_avg = (k_array-1).*NonTF_z_avg./k_array + NonTF_z./k_array;
end

Prop_z_avg  = single(Prop_z_avg);
Dissp_z_avg = single(Dissp_z_avg);
NonTp_z_avg = single(NonTp_z_avg);
ProF_avg  = single(ProF_avg);
DissF_avg = single(DissF_avg);
NonTF_avg = single(NonTF_avg);
ProF_z_avg  = single(ProF_z_avg);
DissF_z_avg = single(DissF_z_avg);
NonTF_z_avg = single(NonTF_z_avg);



savename = strcat('avg_PDN_pF',num2str(Retau),'.mat');
save(savename,'Prop_z_avg','Dissp_z_avg','NonTp_z_avg','ProF_avg','DissF_avg','NonTF_avg','ProF_z_avg','DissF_z_avg','NonTF_z_avg','-v7.3')

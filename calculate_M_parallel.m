function calculate_M_parallel(Retau,k_edge1,k_edge2,read_array,jobid,workers) 
myCluster = parcluster('local');
myCluster.NumWorkers = workers;
parpool(workers)

fft_x = @(x) fft(x, [], 3) / size(x, 3);
fft_y = @(x) fft(x, [], 2) / size(x, 2);
fft_xy = @(x) fft_x(fft_y(x));

load(['full',num2str(Retau),'_mean.mat']); channelRe=channelRe; dkx=dkx; dky=dky; kx_array=kx_array; ky_array=ky_array; Lx=Lx; Ly=Ly; nx=nx; ny=ny; nz=nz; nzDNS=nzDNS; xp=xp; xu=xu; yp=yp; yv=yv; zp=zp; zw=zw;
u_avg_reference = channelRe.Up;

[Diff,zc] = cheb(nz);
Diff = Diff(2:end-1,2:end-1);

if Retau == 180
    dirname = 'grid_180/112x112x150';
else
    dirname = 'grid_590/384x384x500';
end

kx = [0:dkx:(k_edge1*dkx),-(k_edge1*dkx):dkx:-dkx];
ky = [0:dky:(k_edge1*dky)];
kx_array = repelem(kx,length(ky));
ky_array = repmat(ky,1,length(kx));
kx_array = kx_array.';
ky_array = ky_array.';

kx_0posi = [0:dkx:k_edge2*dkx];
ky_0posi = [0:dky:k_edge2*dky];
kx_posi  = kx_0posi(2:end);
ky_posi  = ky_0posi(2:end);
N_wavenumber = (length(kx_0posi)*length(ky_0posi)-1)*(length(kx_0posi)*length(ky_0posi)-2)/2;
M_4d_avg     = zeros(length(ky_0posi),length(kx_0posi),length(ky_0posi),length(kx_0posi));
% M_5d_avg     = zeros(nz-1,length(ky_0posi),length(kx_0posi),length(ky_0posi),length(kx_0posi));
[kx_m,ky_m] = getkm(kx_array,ky_array,nx,ny,dkx,dky);

for k_array = 1:length(read_array)
    k_array
    u = zeros(nzDNS+2,ny,nx);
    v = zeros(nzDNS+2,ny,nx);
    w = zeros(nzDNS+1,ny,nx);
    loadname = strcat('../ChanFast/',dirname,'/outputdir/u_it',num2str(read_array(k_array),'%.0f'),'.dat');
    u(2:end-1,:,:) = read_bin(loadname, [nzDNS, ny, nx]);
    u = interp1(zp,u,zc) - u_avg_reference;
    u = u(2:end-1,:,:);
    u = permute(u,[3 2 1]);
    u = interp1(xu,u,xp,'linear','extrap');
    u = permute(u,[3 2 1]);
    loadname = strcat('../ChanFast/',dirname,'/outputdir/v_it',num2str(read_array(k_array),'%.0f'),'.dat');
    v(2:end-1,:,:) = read_bin(loadname, [nzDNS, ny, nx]);
    v = interp1(zp,v,zc);
    v = v(2:end-1,:,:);
    v = permute(v,[2 1 3]);
    v = interp1(yv,v,yp,'linear','extrap');
    v = permute(v, [2 1 3]);
    loadname = strcat('../ChanFast/',dirname,'/outputdir/w_it',num2str(read_array(k_array),'%.0f'),'.dat');
    w(2:end,:,:) = read_bin(loadname, [nzDNS, ny, nx]);
    w = interp1(zw,w,zc);
    w = w(2:end-1,:,:);
    
    [~,WEIGHT] = clenCurt(nz);
    WEIGHT = WEIGHT(2:nz);
    
    u_F  = fft_xy(u);
    v_F  = fft_xy(v);
    w_F  = fft_xy(w);
    
    u_F(:,1,1) = 0;
    v_F(:,1,1) = 0;
    w_F(:,1,1) = 0;
    
    [duF_dx,duF_dy,duF_dz]    = get_3d(u_F,Diff,kx_m,ky_m);
    [dvF_dx,dvF_dy,dvF_dz]    = get_3d(v_F,Diff,kx_m,ky_m);
    [dwF_dx,dwF_dy,dwF_dz]    = get_3d(w_F,Diff,kx_m,ky_m); 

    [M_4d] = calculate_M_abkxky(kx_0posi,ky_0posi,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky,WEIGHT);
    
    M_4d_avg     = (k_array-1).*M_4d_avg./k_array + M_4d./k_array;
    % M_5d_avg     = (k_array-1).*M_5d_avg./k_array + M_5d./k_array;
end

delete(gcp('nocreate'))

% M_5d_avg     = single(M_5d_avg);
M_4d_avg     = single(M_4d_avg);
kx_0posi     = single(kx_0posi);
ky_0posi     = single(ky_0posi);

savename = strcat('M_',num2str(Retau),'_',num2str(jobid),'.mat');
save(savename,'M_4d_avg','kx_0posi','ky_0posi','-v7.3')

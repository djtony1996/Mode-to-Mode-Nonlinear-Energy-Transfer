clear,clc

Retau = 590;
load(strcat('full',num2str(Retau),'_mean.mat'));
U    = channelRe.Up;
dUdz = channelRe.Up_diff1(2:end-1);

load(strcat('avg_p_physical_space',num2str(Retau),'.mat'));

if Retau == 180
    dirname = 'grid_180/112x112x150';
    kwave_max = 50;
else
    dirname = 'grid_590/384x384x500';
    kwave_max = 150;
end

[Diff,zc] = cheb(nz);
Diff = Diff(2:end-1,2:end-1);
[~,WEIGHT] = clenCurt(nz);
WEIGHT = WEIGHT(2:end-1);

kx = [0:dkx:(kwave_max*dkx),-(kwave_max*dkx):dkx:-dkx];
ky = [0:dky:(kwave_max*dky)];
kx_array = repelem(kx,length(ky));
ky_array = repmat(ky,1,length(kx));
kx_array = kx_array.';
ky_array = ky_array.';

read_array = [100000:10000:890000];

Prop_avg    = zeros(nz-1,1);
Dissp_avg   = zeros(nz-1,1);
NonTp_avg   = zeros(nz-1,1);
P_diffp_avg = zeros(nz-1,1);
V_diffp_avg = zeros(nz-1,1);

for k_array = 1:length(read_array)
    k_array
    u = zeros(nzDNS+2,ny,nx);
    v = zeros(nzDNS+2,ny,nx);
    w = zeros(nzDNS+1,ny,nx);
    loadname = strcat('../ChanFast/',dirname,'/outputdir/u_it',num2str(read_array(k_array),'%.0f'),'.dat');
    u(2:end-1,:,:) = read_bin(loadname, [nzDNS, ny, nx]);
    u = interp1(zp,u,zc);
    u = u - U;
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
    loadname = strcat('../ChanFast/',dirname,'/outputdir/p_it',num2str(read_array(k_array),'%.0f'),'.dat');
    p = read_bin(loadname, [nzDNS, ny, nx]);
    p = interp1(zp(2:end-1),p,zc(2:end-1)) - p_avg;

    [Prop,Dissp,NonTp,P_diffp,V_diffp] = get_five_energy_physicalspace(u,v,w,p,kx_array,ky_array,nx,ny,dkx,dky,Diff,dUdz,Retau);
    Prop    = mean(Prop,[2 3]);
    Dissp   = mean(Dissp,[2 3]);
    NonTp   = mean(NonTp,[2 3]);
    P_diffp = mean(P_diffp,[2 3]);
    V_diffp = mean(V_diffp,[2 3]);

    % calculate the average quantities
    Prop_avg    = (k_array-1).*Prop_avg./k_array + Prop./k_array;
    Dissp_avg   = (k_array-1).*Dissp_avg./k_array + Dissp./k_array;
    NonTp_avg   = (k_array-1).*NonTp_avg./k_array + NonTp./k_array;
    P_diffp_avg = (k_array-1).*P_diffp_avg./k_array + P_diffp./k_array;
    V_diffp_avg = (k_array-1).*V_diffp_avg./k_array + V_diffp./k_array;
end


savename = strcat('avg_energy_physical_space',num2str(Retau),'.mat');
save(savename,'Prop_avg','Dissp_avg','NonTp_avg','P_diffp_avg','V_diffp_avg','-v7.3')

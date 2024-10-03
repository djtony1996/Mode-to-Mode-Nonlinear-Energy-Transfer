function [Prop,Dissp,NonTp,P_diffp,V_diffp] = get_five_energy_physicalspace(u,v,w,p,kx_array,ky_array,nx,ny,dkx,dky,Diff,dUdz,Retau)

fft_x = @(x) fft(x, [], 3) / size(x, 3);
fft_y = @(x) fft(x, [], 2) / size(x, 2);
fft_xy = @(x) fft_x(fft_y(x));

ifft_x = @(x) ifft(x, [], 3) * size(x, 3);
ifft_y = @(x) ifft(x, [], 2) * size(x, 2);
ifft_xy = @(x) ifft_x(ifft_y(x));

uu = u.*u;
vv = v.*v;
ww = w.*w;
up = u.*p;
vp = v.*p;
wp = w.*p;

u_F  = fft_xy(u);
v_F  = fft_xy(v);
w_F  = fft_xy(w);
uu_F = fft_xy(uu);
vv_F = fft_xy(vv);
ww_F = fft_xy(ww);
up_F = fft_xy(up);
vp_F = fft_xy(vp);
wp_F = fft_xy(wp);

[kx_m,ky_m] = getkm(kx_array,ky_array,nx,ny,dkx,dky);
[duF_dx,duF_dy,duF_dz]    = get_3d(u_F,Diff,kx_m,ky_m);
[dvF_dx,dvF_dy,dvF_dz]    = get_3d(v_F,Diff,kx_m,ky_m);
[dwF_dx,dwF_dy,dwF_dz]    = get_3d(w_F,Diff,kx_m,ky_m);
% [dpF_dx,dpF_dy,dpF_dz]    = get_3d(p_F,Diff,kx_m,ky_m);
[dupF_dx,~,~]    = get_3d(up_F,Diff,kx_m,ky_m);
[~,dvpF_dy,~]    = get_3d(vp_F,Diff,kx_m,ky_m);
[~,~,dwpF_dz]    = get_3d(wp_F,Diff,kx_m,ky_m);
[d2uuF_dx2,d2uuF_dy2,d2uuF_dz2]           = get_3d_2nd(uu_F,Diff,kx_m,ky_m);
[d2vvF_dx2,d2vvF_dy2,d2vvF_dz2]           = get_3d_2nd(vv_F,Diff,kx_m,ky_m);
[d2wwF_dx2,d2wwF_dy2,d2wwF_dz2]           = get_3d_2nd(ww_F,Diff,kx_m,ky_m);

Prop = -u.*w.*dUdz;

Dissxp = (-ifft_xy(duF_dx).*ifft_xy(duF_dx) - ifft_xy(duF_dy).*ifft_xy(duF_dy) - ifft_xy(duF_dz).*ifft_xy(duF_dz))./Retau;
Dissyp = (-ifft_xy(dvF_dx).*ifft_xy(dvF_dx) - ifft_xy(dvF_dy).*ifft_xy(dvF_dy) - ifft_xy(dvF_dz).*ifft_xy(dvF_dz))./Retau;
Disszp = (-ifft_xy(dwF_dx).*ifft_xy(dwF_dx) - ifft_xy(dwF_dy).*ifft_xy(dwF_dy) - ifft_xy(dwF_dz).*ifft_xy(dwF_dz))./Retau;
Dissp  = Dissxp + Dissyp + Disszp;

NonTxp = -u.*u.*ifft_xy(duF_dx) - u.*v.*ifft_xy(duF_dy) - u.*w.*ifft_xy(duF_dz);
NonTyp = -v.*u.*ifft_xy(dvF_dx) - v.*v.*ifft_xy(dvF_dy) - v.*w.*ifft_xy(dvF_dz);
NonTzp = -w.*u.*ifft_xy(dwF_dx) - w.*v.*ifft_xy(dwF_dy) - w.*w.*ifft_xy(dwF_dz);
NonTp  = NonTxp + NonTyp + NonTzp;

P_diffp = - ifft_xy(dupF_dx) - ifft_xy(dvpF_dy) - ifft_xy(dwpF_dz);

V_diffxp = ifft_xy(d2uuF_dx2)./2./Retau + ifft_xy(d2uuF_dy2)./2./Retau + ifft_xy(d2uuF_dz2)./2./Retau;
V_diffyp = ifft_xy(d2vvF_dx2)./2./Retau + ifft_xy(d2vvF_dy2)./2./Retau + ifft_xy(d2vvF_dz2)./2./Retau;
V_diffzp = ifft_xy(d2wwF_dx2)./2./Retau + ifft_xy(d2wwF_dy2)./2./Retau + ifft_xy(d2wwF_dz2)./2./Retau;
V_diffp  = V_diffxp + V_diffyp + V_diffzp;

end



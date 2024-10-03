% calculate the three energy: production (Pro), Diss (dissipation), NonT
% (nonlinear transfer) in Fourier space
function [Pro, Diss, NonT, Pro_3, Diss_3, NonT_3] = get_three_energy(u,v,w,kx_array,ky_array,nx,ny,nz,dkx,dky,Retau,dUdz,Diff)

[kx_m,ky_m] = getkm(kx_array,ky_array,nx,ny,dkx,dky);

fft_x = @(x) fft(x, [], 3) / size(x, 3);
fft_y = @(x) fft(x, [], 2) / size(x, 2);
fft_xy = @(x) fft_x(fft_y(x));
ifft_x = @(x) ifft(x, [], 3) * size(x, 3);
ifft_y = @(x) ifft(x, [], 2) * size(x, 2);
ifft_xy = @(x) ifft_x(ifft_y(x));    

uu = u .* u;
vv = v .* v;
ww = w .* w;
uv = u .* v;
uw = u .* w;
vw = v .* w;

u_F  = fft_xy(u);
v_F  = fft_xy(v);
w_F  = fft_xy(w);
uu_F = fft_xy(uu);
vv_F = fft_xy(vv);
ww_F = fft_xy(ww);
uv_F = fft_xy(uv);
uw_F = fft_xy(uw);
vw_F = fft_xy(vw);

[duF_dx,duF_dy,duF_dz]    = get_3d(u_F,Diff,kx_m,ky_m);
[dvF_dx,dvF_dy,dvF_dz]    = get_3d(v_F,Diff,kx_m,ky_m);
[dwF_dx,dwF_dy,dwF_dz]    = get_3d(w_F,Diff,kx_m,ky_m);
[duuF_dx,~,~]             = get_3d(uu_F,Diff,kx_m,ky_m);
[~,dvvF_dy,~]             = get_3d(vv_F,Diff,kx_m,ky_m);
[~,~,dwwF_dz]             = get_3d(ww_F,Diff,kx_m,ky_m);
[duvF_dx,duvF_dy,~]       = get_3d(uv_F,Diff,kx_m,ky_m);
[duwF_dx,~,duwF_dz]       = get_3d(uw_F,Diff,kx_m,ky_m);
[~,dvwF_dy,dvwF_dz]       = get_3d(vw_F,Diff,kx_m,ky_m);

Pro_3  = -(dUdz.*u_F) .* conj(w_F) - conj(dUdz.*u_F) .* w_F;
Diss_3 = 2.* (-duF_dx.*conj(duF_dx)-duF_dy.*conj(duF_dy)-duF_dz.*conj(duF_dz)-dvF_dx.*conj(dvF_dx)-dvF_dy.*conj(dvF_dy)-dvF_dz.*conj(dvF_dz)-dwF_dx.*conj(dwF_dx)-dwF_dy.*conj(dwF_dy)-dwF_dz.*conj(dwF_dz))./Retau;
NonT_3 = -u_F.*conj(duuF_dx)-u_F.*conj(duvF_dy)-u_F.*conj(duwF_dz)-v_F.*conj(duvF_dx)-v_F.*conj(dvvF_dy)-v_F.*conj(dvwF_dz)-w_F.*conj(duwF_dx)-w_F.*conj(dvwF_dy)-w_F.*conj(dwwF_dz) - conj(u_F).*duuF_dx-conj(u_F).*duvF_dy-conj(u_F).*duwF_dz-conj(v_F).*duvF_dx-conj(v_F).*dvvF_dy-conj(v_F).*dvwF_dz-conj(w_F).*duwF_dx-conj(w_F).*dvwF_dy-conj(w_F).*dwwF_dz;

% du_dx = ifft_xy(duF_dx); du_dy = ifft_xy(duF_dy); du_dz = ifft_xy(duF_dz);
% dv_dx = ifft_xy(dvF_dx); dv_dy = ifft_xy(dvF_dy); dv_dz = ifft_xy(dvF_dz);
% dw_dx = ifft_xy(dwF_dx); dw_dy = ifft_xy(dwF_dy); dw_dz = ifft_xy(dwF_dz);
% 
% uduF_dx = fft_xy(u.*du_dx); vduF_dy = fft_xy(v.*du_dy); wduF_dz = fft_xy(w.*du_dz);
% udvF_dx = fft_xy(u.*dv_dx); vdvF_dy = fft_xy(v.*dv_dy); wdvF_dz = fft_xy(w.*dv_dz);
% udwF_dx = fft_xy(u.*dw_dx); vdwF_dy = fft_xy(v.*dw_dy); wdwF_dz = fft_xy(w.*dw_dz);
% 
% NonT_3 = -u_F.*conj(uduF_dx) - u_F.*conj(vduF_dy) - u_F.*conj(wduF_dz) - v_F.*conj(udvF_dx) - v_F.*conj(vdvF_dy) - v_F.*conj(wdvF_dz) - w_F.*conj(udwF_dx) - w_F.*conj(vdwF_dy) - w_F.*conj(wdwF_dz) - conj(u_F).*uduF_dx - conj(u_F).*vduF_dy - conj(u_F).*wduF_dz - conj(v_F).*udvF_dx - conj(v_F).*vdvF_dy - conj(v_F).*wdvF_dz - conj(w_F).*udwF_dx - conj(w_F).*vdwF_dy - conj(w_F).*wdwF_dz;

[~,W] = clenCurt(nz);
W = W(2:nz);

Diss = squeeze(pagemtimes(W,Diss_3))./2;
Pro  = squeeze(pagemtimes(W,Pro_3))./2;
NonT = squeeze(pagemtimes(W,NonT_3))./2;

end

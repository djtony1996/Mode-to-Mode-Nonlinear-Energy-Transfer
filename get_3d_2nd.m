% get the 2nd order derivatives in x,y,z directions (in Fourier space x-y)
function [d2fdx2,d2fdy2,d2fdz2] = get_3d_2nd(f,Diff,kx_m,ky_m)
    f = permute(f,[2 3 1]);
    d2fdx2 = - (kx_m.^2) .* f;
    d2fdy2 = - (ky_m.^2) .* f;
    
    f    = permute(f,[3 1 2]);
    d2fdx2 = permute(d2fdx2,[3 1 2]);
    d2fdy2 = permute(d2fdy2,[3 1 2]);
    
    d2fdz2  = pagemtimes(Diff*Diff,f);
end

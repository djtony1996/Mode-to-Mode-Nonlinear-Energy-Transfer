% get the derivatives in x,y,z directions (in Fourier space x-y)
function [dfdx,dfdy,dfdz] = get_3d(f,Diff,kx_m,ky_m)
    f = permute(f,[2 3 1]);
    dfdx = 1i .* kx_m .* f;
    dfdy = 1i .* ky_m .* f;
    
    f    = permute(f,[3 1 2]);
    dfdx = permute(dfdx,[3 1 2]);
    dfdy = permute(dfdy,[3 1 2]);
    
    dfdz  = pagemtimes(Diff,f);
end

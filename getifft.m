% This function is used to implement inverse Fourier transform to get the
% velocity in time domain. 
%
% Basically, it accepts the Fourier coefficients and the corresponding
% wavenumbers. Form the Fourier coefficient matrix and then use the 
% Matlab command 'ifft2'. 
%
% One thing to note is that this Fourier coefficient matrix construction is
% only for the streamwise wavenumber array and the spanwise wavenumber
% array of the form: kx = [kx_min:dkx:kx_max] ky = [ky_min:dky:ky_max]. 
% The order of elements in the arrays does not matter, but the elements 
% in these two arrays should correspond to each other.
%
% For more information about constructing the Fourier coefficient matrix,
% please look at the comments in code 'inverseF1.m' in folder 'channel
% flow/ifft'.
function [velpro, Fcoeffmatrix] = getifft(vel,alpha_array,beta_array,nx,ny,dkx,dky)
    Fcoeffmatrix = zeros(ny,nx);
    velpro = zeros(ny,nx);        % velocity to be calculated using 'ifft2'

    for k_array = 1: length(alpha_array)

        index_kx = round(alpha_array(k_array)/dkx) + 1;
        index_ky = round(beta_array(k_array)/dky) + 1;

        if alpha_array(k_array)<0
            Fcoeffmatrix(index_ky,size(velpro,2)+index_kx) = vel(k_array);
        else
            Fcoeffmatrix(index_ky,index_kx) = vel(k_array);
        end
    end

    Fcoeffmatrix(1,(end:-1:(end-(nx/2-2)))) = conj(Fcoeffmatrix(1,(2:(nx/2)))); % (kx 0)   <->  (-kx 0)
    Fcoeffmatrix((end:-1:(end-(ny/2-2))),1) = conj(Fcoeffmatrix((2:(ny/2)),1)); % (0 ky)   <->  (0 -ky)

    Fcoeffmatrix((end:-1:(end-(ny/2-2))),(end:-1:(end-(nx/2-2)))) = conj(Fcoeffmatrix(2:(ny/2),2:(nx/2))); % (kx ky)  <->  (-kx -ky)
    Fcoeffmatrix(end:-1:(end-(ny/2-2)),2:(nx/2)) = conj(Fcoeffmatrix(2:(ny/2),end:-1:(end-(nx/2-2))));     % (kx -ky) <->  (-kx ky)

    velpro(:,:) = ifft2(Fcoeffmatrix).*size(Fcoeffmatrix,1).*size(Fcoeffmatrix,2);    % pay attention to normalisation

end

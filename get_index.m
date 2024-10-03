% This function is to get the x-coordinates and y-coordinates of the
% wavenumber in 2D Fourier space (MATLAB fft2)
function [index_x, index_y] = get_index(wave_x,wave_y,nx,ny,dkx,dky)

        index_kx = round(wave_x/dkx) + 1;
        index_ky = round(wave_y/dky) + 1;

        if wave_x<0
            index_x = nx + index_kx;
        else
            index_x = index_kx;
        end

        if wave_y<0
            index_y = ny + index_ky;
        else
            index_y = index_ky;
        end
        
end
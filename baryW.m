function b = baryW(N)

b = [1/2; ones(N-1,1); 1/2] .* (-1).^((0:N)');
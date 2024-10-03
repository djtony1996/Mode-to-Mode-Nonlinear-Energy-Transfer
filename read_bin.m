function A = read_bin(filename, dim)
% A = read_bin(filename, [nz, ny ,...]) will read a binary file of this
% size.

fid = fopen(filename);
A = reshape(fread(fid,prod(dim), 'double', 'l'), dim);
fclose(fid);

end
function out = ifftnc(in)
out = in;
for n = 1:ndims(in)
    out = ifftc(out,n);
end
end
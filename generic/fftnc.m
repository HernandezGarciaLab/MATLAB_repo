function out = fftnc(in)
out = in;
for n = 1:ndims(in)
    out = fftc(out,n);
end
end
function out = ifft3c(in)
out = ifftc(ifftc(ifftc(in,1),2),3);
end
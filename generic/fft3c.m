function out = fft3c(in)
out = fftc(fftc(fftc(in,1),2),3);
end
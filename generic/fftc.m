function y = fftc(x,dim)
y = 1/sqrt(size(x,dim))*fftshift(fft(fftshift(x,dim),[],dim),dim);
end
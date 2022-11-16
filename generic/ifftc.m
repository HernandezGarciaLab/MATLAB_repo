function y = ifftc(x,dim)
y = sqrt(size(x,dim))*ifftshift(ifft(ifftshift(x,dim),[],dim),dim);
end
function out = fftc(in,dim)

    % Set default dimensions
    if nargin<2 || isempty(dim)
        dim = 1:ndims(in);
    elseif any(dim(:) > ndims(in)) || any(dim(:) < 1)
        error('dimensions out of range');
    end
    
    % Define fourier transform with scaling and shifts
    fftc1d = @(x,d) 1/sqrt(size(x,d))*fftshift(fft(fftshift(x,d),[],d),d);
    
    % Fourier transform along each requested dimension
    out = in;
    for n = dim
        out = fftc1d(out,n);
    end
    
end
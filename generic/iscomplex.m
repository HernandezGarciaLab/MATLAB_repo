function TF = iscomplex(val)
    TF = any(imag(val),'all') > 0;
end
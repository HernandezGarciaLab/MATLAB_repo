function varargout = quickgrid(fov,dim)

    % Check that fov and dim are same size
    if numel(fov) ~= numel(dim)
        error('fov and dim must be same size');
    end

    % Check that number of grid dimensions in output is equal to number of
    % input dimensions
    if nargout ~= numel(fov)
        error('number of arguments out must be equal to number of dimenstions in fov');
    end

    % Loop through dimensions until each is accounted for in grid
    outcode = 'varargout{1}';
    incode = 'linspace(-fov(1)/2,fov(1)/2,dim(1))';
    varargout = cell(1,nargout);
    if numel(fov) > 1
        for n = 2:numel(fov)
            outcode = [outcode, sprintf(',varargout{%d}',n)];
            incode = [incode, sprintf(',linspace(-fov(%d)/2,fov(%d)/2,dim(%d))',n,n,n)];
        end
    end
    eval(sprintf('[%s] = ndgrid(%s);',outcode,incode));

end


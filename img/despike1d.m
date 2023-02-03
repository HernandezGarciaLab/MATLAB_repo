function data_corr = despike1d(data,dim,varargin)

if nargin < 2 || isempty(dim)
    dim = 1;
end

defaults = struct( ...
    'linked',           {1:size(data,dim)}, ...
    'thresh',           3, ...
    'outliermethod',    'median' ...
    );

args = vararginparser(defaults,varargin{:});

% Determine permutation by swapping first dimension with desired
permutation = 1:ndims(data);
permutation(dim) = 1;
permutation(1) = dim;

% And permute/reshape the data
data = permute(data,permutation);
datashape = size(data);
data = reshape(data,size(data,1),[]);
data_corr = data;

% Loop through sets of indicies
nbad = 0;
if ~iscell(args.linked)
    args.linked = {args.linked};
end
for s = 1:length(args.linked)
    data_s = data(args.linked{s},:);
    [badx,bady] = find(isoutlier(abs(data_s),args.outliermethod,1));
    nbad = nbad + length(badx);

    % Loop through all unique points that contain at least 1 outlier
    for y = reshape(unique(bady),1,[])

        % Determine & seperate bad/good indicies for point
        x1 = badx(bady == y);
        x0 = 1:length(args.linked{s});
        x0(x1) = [];

        % Interpolate using least squares regression
        data_s(x1,y) = x1(:).^[1,0] * pinv(x0(:).^[1,0]) * data_s(x0,y);

    end

    % Apply 
    data_corr(args.linked{s},:) = data_s;

end

% Fix permutation and reshape
data_corr = reshape(data_corr,datashape);
data_corr = permute(data_corr,permutation);

fprintf('despike1d removed %d points\n',nbad);

end


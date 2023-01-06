function [Sk,Si] = kspace2dphantom(k,x,varargin)
% function [Sk,Si] = kspace2dphantom(k,x,varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to analytically derive kspace signal at k-space
%   sampling points of digital 2d phantoms using method described in Guan 
%   Koay, C., Sarlls, J.E., Özarslan, E., (2007) Three-dimensional
%   analytical magnetic resonance imaging phantom in the Fourier domain,
%   Magn. Reson. Med. 58(2):430-436, https://doi.org/10.1002/mrm.21292
%
%
% Notes:
%   - default values in help message may not be up to date - check defaults
%       structure under the function header
%
% Dependencies:
%   - matlab default path
%       - can be restored by typing 'restoredefaultpath'
%   - umasl
%       - github: fmrifrey/umasl
%       - umasl/matlab/ and subdirectories must be in current path
%
% Static input arguments:
%   - k:
%       - kspace sampling points
%       - Nx2 float/double matrix containing [kx,ky] sampling points
%       - default is empty
%   - x:
%       - image space sampling points, for comparing reconstruction to
%           truth
%       - Nx2 float/double matrix containing [x,y] sampling points
%       - default is empty
%
% Variable input arguments (type 'help varargin' for usage info):
%   - 'e':
%       - image ellipsoid matrix
%       - Can be an Nx6 matrix containing custom ellipsoid info, where the
%           columns are: [intensity, x size, y size, x offset, y offset,
%           z rotation (deg)] respectively, where dimensions are fractions
%           of fov
%       - Can be a string name of a preset 2d phantom: 'shepp-logan' or
%           'modified-shepp-logan'
%       - default is 'modified-shepp-logan'
%   - 'size':
%       - field of view of phantom image (cm)
%       - 1x2 double/float array describing phantom's fov
%       - default is [24,24]
%
% Function output:
%   - Sk:
%       - kspace data
%       - Nx1 vector containing data at Nx2 kspace sampling points (k)
%   - Si:
%       - image space data
%       - Nx1 vector containing data at Nx2 image space sampling points (x)
%

    % Define default arguments
    if nargin < 1 || isempty(k)
        k = [];
    end

    if nargin < 2 || isempty(x)
        x = [];
    end
    
    defaults = struct( ...
        'e', 'modified-shepp-logan', ...
        'size', [24 24] ...
        );
    
    % Parse through variable inputs using matlab's built-in input parser
    args = vararginparser(defaults, varargin{:});

    % Get ellipsoid parameter matrix
    if ischar(args.e)
        % If e is a phantom name
        args.e = lower(args.e);
        switch args.e
            case 'shepp-logan'
                args.e = shepp_logan;
            case 'modified-shepp-logan'
                args.e = modified_shepp_logan;
            otherwise
                error('invalid phantom name: %s', args.e);
        end
    elseif size(args.e, 2) ~= 6
        % If e is not an Nx6 matrix
        error('e must either be a valid phantom name or an Nx10 matrix');
    end
    
    Sk = zeros(size(k,1),1);
    Si = zeros(size(x,1),1);
    
    for ei = 1:size(args.e,1)
        
        p = args.e(ei,1); % ellipsoid amplitude
        ab = args.size(:)'/2.*args.e(ei,2:3); % ellipsoid radii
        d = args.size(:)'/2.*args.e(ei,4:5); % ellipsoid central displacement
        psi = args.e(ei,6)*pi/180; % z Euler angle
       
        % Euler rotation matrix
        cpsi = cos(psi);
        spsi = sin(psi);
        A = [cpsi, spsi; -spsi, cpsi];
       
        % Calculate kspace signal
        if size(k,2) == 2
            kn0 = k(vecnorm(k,2,2)~=0,:);
            ktd = kn0*d';
            K = vecnorm(kn0*A'.*ab,2,2);
            Sk(vecnorm(k,2,2)~=0) = Sk(vecnorm(k,2,2)~=0) + ...
                p * abs(det(A')) * prod(ab) * exp(-1i * 2*pi*ktd) .* ...
                besselj(1, 2 * pi * K) ./ (pi * K);
            
            k0 = k(vecnorm(k,2,2)==0,:);
            ktd = k0*d';
            K = vecnorm(k0*A'.*ab,2,2);
            Sk(vecnorm(k,2,2)==0) = Sk(vecnorm(k,2,2)==0) + ...
                p * abs(det(A')) * prod(ab) * exp(-1i * 2*pi*ktd) .* ...
                (1 - pi^2/2*K.^2 + pi^4/12*K.^4); 
        elseif ~isempty(k)
            error('k must be an Nx2 matrix');
        end
        
        % Calculate image space signal
        if size(x,2) == 2
            Si( sum(((x*A' - d) ./ ab).^2, 2) <= 1 ) = ...
                Si( sum(((x*A' - d) ./ ab).^2, 2) <= 1 ) + p;
        elseif ~isempty(x)
            error('x must be an Nx2 matrix');
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Default head phantoms:   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function e = shepp_logan

e = modified_shepp_logan;
e(:,1) = [1 -.98 -.02 -.02 .01 .01 .01 .01 .01 .01];

end
      
function e = modified_shepp_logan
%
%   This head phantom is the same as the Shepp-Logan except 
%   the intensities are changed to yield higher contrast in
%   the image.  Taken from Toft, 199-200.
%      
%         A    a     b    x0    y0    phi
%        ---------------------------------
e = [     1   .69   .92    0     0     0   
        -.8  .6624 .8740   0  -.0184   0
        -.2  .1100 .3100  .22    0    -18
        -.2  .1600 .4100 -.22    0     18
         .1  .2100 .2500   0    .35    0
         .1  .0460 .0460   0    .1     0
         .1  .0460 .0460   0   -.1     0
         .1  .0460 .0230 -.08  -.605   0 
         .1  .0230 .0230   0   -.606   0
         .1  .0230 .0460  .06  -.605   0   ];
       
end

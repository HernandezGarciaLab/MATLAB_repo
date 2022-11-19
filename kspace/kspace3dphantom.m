function [Sk,Si] = kspace3dphantom(k,x,varargin)
% function [Sk,Si] = kspace3dphantom(k,x,varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to analytically derive kspace signal at k-space
%   sampling points of digital 3d phantoms using method described in Guan 
%   Koay, C., Sarlls, J.E., Özarslan, E., (2007) Three-dimensional
%   analytical magnetic resonance imaging phantom in the Fourier domain,
%   Magn. Reson. Med. 58:2:430-436, https://doi.org/10.1002/mrm.21292
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
%       - Nx3 float/double matrix containing [kx,ky,kz] sampling points
%       - default is empty
%   - x:
%       - image space sampling points, for comparing reconstruction to
%           truth
%       - Nx3 float/double matrix containing [x,y,z] sampling points
%       - default is empty
%
% Variable input arguments (type 'help varargin' for usage info):
%   - 'e':
%       - image ellipsoid matrix
%       - Can be an Nx10 matrix containing custom ellipsoid info, where the
%           columns are: [intensity, x size, y size, c size, x offset, y
%           offset, z offset, x rotation (deg), y rotation (deg), z rotatin
%           (deg)] respectively, where dimensions are fractions of fov
%       - Can be a string name of a preset 3d phantom: 'shepp-logan',
%           'modified-shepp-logan', or 'ye-yu-wang'
%       - default is 'modified-shepp-logan'
%   - 'size':
%       - field of view of phantom image (cm)
%       - 1x3 double/float array describing phantom's fov
%       - default is [24,24,24]
%
% Function output:
%   - Sk:
%       - kspace data
%       - Nx1 vector containing data at Nx3 kspace sampling points (k)
%   - Si:
%       - image space data
%       - Nx1 vector containing data at Nx3 image space sampling points (x)
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
        'size', [24 24 24] ...
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
            case 'yu-ye-wang'
                args.e = yu_ye_wang;
            otherwise
                error('invalid phantom name: %s', args.e);
        end
    elseif size(args.e, 2) ~= 10
        % If e is not an Nx10 matrix
        error('e must either be a valid phantom name or an Nx10 matrix');
    end
    
    Sk = zeros(size(k,1),1);
    Si = zeros(size(x,1),1);
    
    for ei = 1:size(args.e,1)
        
        p = args.e(ei,1); % ellipsoid amplitude
        abc = args.size(:)'/2.*args.e(ei,2:4); % ellipsoid radii
        d = args.size(:)'/2.*args.e(ei,5:7); % ellipsoid central displacement
        phi = args.e(ei,8)*pi/180; % x Euler angle
        theta = args.e(ei,9)*pi/180; % y Euler angle
        psi = args.e(ei,10)*pi/180; % z Euler angle
       
        % Euler rotation matrix
        cphi = cos(phi);
        sphi = sin(phi);
        ctheta = cos(theta);
        stheta = sin(theta);
        cpsi = cos(psi);
        spsi = sin(psi);
        A = [cpsi*cphi-ctheta*sphi*spsi   cpsi*sphi+ctheta*cphi*spsi  spsi*stheta;
            -spsi*cphi-ctheta*sphi*cpsi  -spsi*sphi+ctheta*cphi*cpsi cpsi*stheta;
            stheta*sphi                  -stheta*cphi                ctheta];
       
        % Calculate kspace signal
        if size(k,2) == 3
            ktd = -k*d';
            K = vecnorm(k*A'.*abc,2,2);
            Sk = Sk + ...
                p * abs(det(A')) * prod(abc) * exp(-1i * 2*pi*ktd) .* ...
                ( sin(2*pi*K) - 2*pi*K.*cos(2*pi*K) ) ./ (2*pi^2*K.^3);
        elseif ~isempty(k)
            error('k must be an Nx3 matrix');
        end
        
        % Calculate image space signal
        if size(x,2) == 3
            Si( sum(((x*A' - d) ./ abc).^2, 2) <= 1 ) = ...
                Si( sum(((x*A' - d) ./ abc).^2, 2) <= 1 ) + p;
        elseif ~isempty(x)
            error('x must be an Nx3 matrix');
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
%         A      a     b     c     x0      y0      z0    phi  theta    psi
%        -----------------------------------------------------------------
e =    [  1  .6900  .920  .810      0       0       0      0      0      0
        -.8  .6624  .874  .780      0  -.0184       0      0      0      0
        -.2  .1100  .310  .220    .22       0       0    -18      0     10
        -.2  .1600  .410  .280   -.22       0       0     18      0     10
         .1  .2100  .250  .410      0     .35    -.15      0      0      0
         .1  .0460  .046  .050      0      .1     .25      0      0      0
         .1  .0460  .046  .050      0     -.1     .25      0      0      0
         .1  .0460  .023  .050   -.08   -.605       0      0      0      0
         .1  .0230  .023  .020      0   -.606       0      0      0      0
         .1  .0230  .046  .020    .06   -.605       0      0      0      0 ];
       
end       

function e = yu_ye_wang
%
%   Yu H, Ye Y, Wang G, Katsevich-Type Algorithms for Variable Radius Spiral Cone-Beam CT
%      
%         A      a     b     c     x0      y0      z0    phi  theta    psi
%        -----------------------------------------------------------------
e =    [  1  .6900  .920  .900      0       0       0      0      0      0
        -.8  .6624  .874  .880      0       0       0      0      0      0
        -.2  .4100  .160  .210   -.22       0    -.25    108      0      0
        -.2  .3100  .110  .220    .22       0    -.25     72      0      0
         .2  .2100  .250  .500      0     .35    -.25      0      0      0
         .2  .0460  .046  .046      0      .1    -.25      0      0      0
         .1  .0460  .023  .020   -.08    -.65    -.25      0      0      0
         .1  .0460  .023  .020    .06    -.65    -.25     90      0      0
         .2  .0560  .040  .100    .06   -.105    .625     90      0      0
        -.2  .0560  .056  .100      0    .100    .625      0      0      0 ];
       
end
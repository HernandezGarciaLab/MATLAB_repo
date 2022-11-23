function Wi = pipedcf(G,itrmax)
% function Wi = pipedcf(G,itrmax)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to generate kspace density weights for Gmri
%   reconstruction using Pipe & Menon method, as described in Pipe, J.G.,
%   Menon, P., (1999) Sampling density compensation in MRI: rationale and
%   an iterative numerical solution, Magn. Reson. Med. 41(1):179-86,
%   https://doi.org/10.1002/(sici)1522-2594(199901)41:1%3C179::aid-mrm25%3E3.0.co;2-v
%
%
% Path dependencies:
%   - matlab default path
%       - can be restored by typing 'restoredefaultpath'
%   - umasl
%       - github: fmrifrey/umasl
%       - umasl/matlab/ and subdirectories must be in current path
%   - mirt (matlab version)
%       - github: JeffFessler/mirt
%       - mirt setup must have successfully ran
%
% Static input arguments:
%   - G:
%       - Gmri or Gnufft reconstruction object
%       - no default, required argument
%   - itrmax:
%       - maximum number of iterations to perform
%       - integer describing number of iterations
%       - default is 15
%
% Function output:
%   - Wi:
%       - kspace density weights
%       - Nx1 vector of same length as number of kspace samples
%

    % Set default for itrmax
    if nargin < 2 || isempty(itrmax)
        itrmax = 15;
    end
    
    % If G is a Gmri object, use its Gnufft object
    if isfield(G,'Gnufft')
        G = G.Gnufft;
    end
    
    % Initialize weights to 1 (psf)
    Wi = ones(size(G,1),1);
    
    % Loop through iterations
    for itr = 1:itrmax
        
        % Pipe algorithm: W_{i+1} = W_{i} / (G * (G' * W_{i}))
        d = real( G.arg.st.interp_table(G.arg.st, ...
            G.arg.st.interp_table_adj(G.arg.st, Wi) ) );
        Wi = Wi ./ d;
        
    end
    
    % Normalize weights
    Wi = Wi / sum(abs(Wi));
    
end


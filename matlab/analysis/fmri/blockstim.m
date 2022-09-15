function [x,x_hires] = blockstim(nframes,tstart,toff,ton,tr)
% function [x,x_hires] = blockstim(nframes,tstart,toff,ton,tr)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to create a hemodynamic response regressor to block
%   stimulation fmri experiment
%
%
% Notes:
%   - defaults should be tailored to umasl eprime block motor vis
%       stimulation experiment
%   - default values in help message may not be up to date - check defaults
%       structure under the function header
%
% Dependencies:
%   - matlab default path
%       - can be restored by typing 'restoredefaultpath'
%   - umasl
%       - github: fmrifrey/umasl
%       - umasl/matlab/ and subdirectories must be in current path
%   - spm12
%       - github: spm/spm12
%       - spm12/ and subdirectories must be in current path
%
% Static input arguments:
%   - nframes:
%       - number of temporal frames
%       - integer describing total number of frames
%       - default is 68
%   - tstart:
%       - onset delay for stimulation
%       - double/float describing start time in seconds
%       - default is 0
%   - toff:
%       - duration of block where stimulation is off
%       - double/float describing off time in seconds
%       - default is 30
%   - ton:
%       - duration of block where stimulation is on
%       - double/float describing on time in seconds
%       - default is 30
%   - tr:
%       - temporal frame repetition time
%       - double/float describing TR in seconds
%       - default is 4.5
%
% Function output:
%   - x:
%       - constructed regressor for hemodynamic stimulation response
%       - double/float column vector of nframes x 1
%   - x_hires:
%       - high resolution response waveform for visualization/debugging
%       - double/float column vector of some arbitrary large length
%

    % Determine defaults
    if nargin < 1 || isempty(nframes)
        nframes = 68;
    end
    
    if nargin < 2 || isempty(tstart)
        tstart = 0;
    end
    
    if nargin < 3 || isempty(toff)
        toff = 30;
    end
    
    if nargin < 4 || isempty(ton)
        ton = 30;
    end
    
    if nargin < 5 || isempty(tr)
        tr = 4.5;
    end

    % Determine length of experiment
    tend = (nframes - 1) * tr;
    
    % Construct hires and lores time arrays
    t_hires = 0:1e-3:tend;
    t = 0:tr:tend;
    
    % Construct stimulation waveform
    stim = 1 * (t_hires >= tstart) .* ...
        (mod(t_hires - tstart, ton + toff) > toff);

    % Get upsampled hemodynamic response function using spm
    [~, p] = spm_hrf(tr * 1e-3);
    p(5) = 100;
    [hrf, ~] = spm_hrf(tr * 1e-3, p);
    
    % Convolve hrf with stimulation waveform and correct length
    x_hires = conv(hrf,stim);
    x_hires = x_hires(1:length(t_hires))';
    
    % Get final regressor
    x = interp1(t_hires, x_hires, t);
    
end


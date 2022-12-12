%% Compute the amplitude envelope correlation as defined in Hipp et al. 2012 Nature Neuroscience
% 
%   input:  - virtChan_timeSeries: cell array of EEG time series, each cell
%           should represent a single epoch/trial
%   output: - connMatrix: s x s x e matrix, containing the connectivity
%           matrices for every single epoch/trial. s = number of sources, e = number of epochs
% 
% Cristina Gil, Felix Bott and Stefan Dvoretskii, Technical University of Munich 12.12.2022
% based on brainstorm AEC connectivity
function connMatrix = aecConnectivity(virtChan_timeSeries)
    nEpochs = length(virtChan_timeSeries.trial);
    nVoxel = length(virtChan_timeSeries.label);
    connMatrix = zeros(nVoxel, nVoxel, nEpochs);
    tol = single(10e-9);
    for iEpoch = 1:nEpochs
        tic;
        % Analytic signal of every virtual time series with Hilbert transform
        HA = hilbert(virtChan_timeSeries.trial{iEpoch}')';
        R = zeros(nVoxel, nVoxel);
        for iSeed = 1:nVoxel
            seed = HA(iSeed,:); % complex-valued signal
            % Orthogonalize all the signals w.r.t. the seed
            HAo = imag(HA.*(conj(seed)./abs(seed))); % this is a real-valued matrix (absolute values)
            % Avoid rounding errors
            HAo(abs(HAo./abs(HA))<2*eps)=0;
            % Correlation of log-transformed Power Envelopes
            c = corr(log(abs(seed).^2 + tol)', log(abs(HAo).^2 + tol)');
            R(iSeed, :) = c;
        end
        % Average correlation in both directions (X->Y and Y->X)
        connMatrix(:,:,iEpoch) = (R+R')/2;
        t = toc;
        disp(['Computing AEC of trial ', num2str(iEpoch) ' took ', num2str(t ,'%.2f'), ' seconds'])
    end
    connMatrix = connMatrix ./ 0.577; % normalization because of under-estimation through orthogonalization, see Hipp et al. 2012 Nature Neuroscience
end
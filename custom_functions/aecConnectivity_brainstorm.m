%% Compute the amplitude envelope correlation as defined in Hipp et al. 2012 Nature Neuroscience
% 
%   input:  - virtChan_timeSeries: cell array of EEG time series, each cell
%           should represent a single epoch/trial
%   output: - connMatrix: s x s x e matrix, containing the connectivity
%           matrices for every single epoch/trial. s = number of sources, e = number of epochs
% 
% Cristina Gil and Stefan Dvoretskii, Technical University of Munich 20.04.2022
% based on brainstorm AEC connectivity
function connMatrix = aecConnectivity_brainstorm(virtChan_timeSeries)
    nTrial = length(virtChan_timeSeries.trial);
    nVoxel = length(virtChan_timeSeries.label);
    connMatrix = zeros(nVoxel, nVoxel, nTrial);
    for iTrial = 1:nTrial
        tic;
        % Estimate envelope of the signals with Hilbert transform
        HA = hilbert(virtChan_timeSeries.trial{iTrial}')';
        R = zeros(nVoxel, nVoxel);
        for iSeed = 1:nVoxel
            % Orthogonalize all the signals w.r.t. the seed
            HAo = imag(bsxfun(@times, HA, conj(HA(iSeed,:))./abs(HA(iSeed,:))));
             % avoid rounding errors
            HAo(abs(HAo./abs(HA))<2*eps)=0;
            % Compute correlation coefficients between the seed and the rest of the signals
            R(iSeed, :) = correlate_dims(abs(HA(iSeed,:)), abs(HAo), 2);
        end
        % Average correlation in both directions (X->Y and Y->X)
        connMatrix(:,:,iTrial) = (R+R')/2;
        t = toc;
        disp(['Computing AEC of trial ', num2str(iTrial) ' took ', num2str(t ,'%.2f'), ' seconds'])
    end
    connMatrix = connMatrix ./ 0.577; % normalization because of under-estimation through orthogonalization, see Hipp et al. 2012 Nature Neuroscience
end

function R = correlate_dims(A, B, dim)
    A = bsxfun( @minus, A, mean( A, dim) );
    B = bsxfun( @minus, B, mean( B, dim) );
    A = normr(A);
    B = normr(B);
    R = sum(bsxfun(@times, A, B), dim);
end

function x = normr(x)
    n = sqrt(sum(x.^2,2));
    x(n~=0,:) = bsxfun(@rdivide, x(n~=0,:), n(n~=0));
    x(n==0,:) = 1 ./ sqrt(size(x,2));
end
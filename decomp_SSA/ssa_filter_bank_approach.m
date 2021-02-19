%-------------------------------------------------------------------------------
% ssa_filter_bank_approach: SSA implemented using a filter-bank approach (see [1-3])
%
% Syntax: [y, x_components, l, U] = ssa_filter_bank_approach(x, L, weights, db_plot)
%
% Inputs: 
%     x       - input signal (1D vector, length-N)
%     L       - embedding dimension (labelled 'M' in [1])
%     weights - vector of weigths when reconstructing signal (1D vector, length-L) 
%     db_plot - plot components with eigenvectors (true/false)
%
% Outputs: 
%     y            - reconstructed signal (1D vector, length-N)
%     x_components - all components from SSA (L x N matrix)
%     l            - eigenvalues (1D vector, length-L)
%     U            - eigenvectors (L x L matrix)
% 
% Example:
%     x = randn(1, 1000);
%     ssa_filter_bank_approach(x, 20, [], true);
% 
%
% [1] Harris, T. J., & Yuan, H. (2010). Filtering and frequency interpretations of
% Singular Spectrum Analysis. Physica D: Nonlinear Phenomena, 239(20–22), 1
% 958–1967. http://doi.org/10.1016/j.physd.2010.07.005
%
% [2] Vautard, R., Yiou, P., & Ghil, M. (1992). Singular-spectrum analysis: A toolkit for
% short, noisy chaotic signals. Physica D: Nonlinear Phenomena, 58(1–4),
% 95–126. http://doi.org/10.1016/0167-2789(92)90103-T
% 
% [3] O'Toole JM. Dempsey EM, Boylan GB (2018) 'Extracting transients from cerebral
% oxygenation signals of preterm infants: a new singular-spectrum analysis method' in
% Int Conf IEEE Eng Med Biol Society (EMBC), IEEE, pp. 5882--5885 
% [https://doi.org/10.1109/EMBC.2018.8513523]


% John M. O' Toole, University College Cork
% Started: 04-11-2016
%
% last update: Time-stamp: <2021-02-16 11:08:09 (otoolej)>
%-------------------------------------------------------------------------------
function [y, x_components, l, U] = ssa_filter_bank_approach(x, L, weights, db_plot)
if(nargin<3 || isempty(weights)), weights = []; end
if(nargin < 4 || isempty(db_plot)), db_plot = false; end



x = x(:)';
N = length(x);


%---------------------------------------------------------------------
% 1. lag matrix
%---------------------------------------------------------------------

% set this to true if MEX files is available (will produce ~ x2 speed-up):
use_mex_file = false;
if(use_mex_file)
    % converted to C code as faster:
    x_lag = fast_lag_matrix_mex(x, int32(L));
    
else
    K = N-L+1;
    L_1 = L-1;
    x_lag = NaN(L, K);

    for k = 1:K
        x_lag(:, k) = x( k:(k+L_1) );
    end
    
end




%---------------------------------------------------------------------
% 2. correlation matrix
%---------------------------------------------------------------------
R = x_lag*x_lag';


%---------------------------------------------------------------------
% 3. eigenvalue decomposition
%---------------------------------------------------------------------
[U, l] = eig(R);
l = diag(l);

[l, il] = sort(l, 'descend');
U = U(:, il);


%---------------------------------------------------------------------
% 4. implement zero-phase filtering using eigenvectors as filters
%---------------------------------------------------------------------
x_components = zeros(L, N);

% implementation equals matlabs 'filtfilt' (when USE_FILTFILT) but faster

% extrapolation signal to migate end-effects from filtering:
USE_FILTFILT = 0;
if(USE_FILTFILT)
    % negative mirror extrapolation (filtfilt implementation):
    % 'reflection method'
    Px = fft([2*x(1)-fliplr(x(2:L+1)) x 2*x(end)-fliplr(x(N-L:end-1))]);    
else
    % or just periodically extend the signal backwards:
    Px = fft([fliplr(x(2:L+1)) x fliplr(x(N-L:end-1))]);
end

in = 1:L;
% only do for non-zero weights:
if(~isempty(weights))
    in = find(weights);    
end
H = fft(U(:, in), length(Px)); 
xf = ifft( bsxfun(@times, abs(H).^2, Px.') )./L;
x_components(in, :) = xf((L+1):(N+L), :).';

if(~isempty(weights))
    x_components = bsxfun(@times, x_components, weights');
end
y = sum(x_components);
    


%---------------------------------------------------------------------
% 5. plot
%---------------------------------------------------------------------
if(db_plot)
    plot_ssa_components(x_components, U, 11);
end




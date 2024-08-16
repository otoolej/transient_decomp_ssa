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

%-------------------------------------------------------------------------------
% trim_transient_component: trim around the impulses/transients 
%
% Syntax: y_trim = trim_transient_component(y)
%
% Inputs: 
%     y - 
%
% Outputs: 
%     y_trim - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 22-04-2021
%
% last update: Time-stamp: <2021-04-27 16:46:16 (otoolej)>
%-------------------------------------------------------------------------------
function y_trim = trim_transient_component(y, fs, db_plot)
if(nargin < 2 || isempty(fs)), fs = 1 / 6; end
if(nargin < 3 || isempty(db_plot)), db_plot = false; end

    
%---------------------------------------------------------------------
% 0. parameters
%---------------------------------------------------------------------
% soft-thresholding parameter:
L_soft = 8;
% smoothing collar for mask:
L_win = 10 * 60 * (1/6);


%---------------------------------------------------------------------
% 1. generate mask from soft threshold
%---------------------------------------------------------------------
y_soft = soft_thres(y, L_soft);

% then convert to discrete [-1, 0, 1] sequence:
y_hsoft = y_soft;
y_hsoft(y_hsoft > 0) = 1;
y_hsoft(y_hsoft < 0) = -1;

% TEST: exclude up-ward transients:
y_hsoft(y_hsoft > 0) = 0;

% find edges and expand with smoothing function:
mask = expand_mask_edges(y_hsoft, 100);

%---------------------------------------------------------------------
% 2. trim by masking signal
%---------------------------------------------------------------------
y_trim = y .* mask;



if(db_plot)
    set_figure(3); 
    plot(y);
    plot(y_soft);
    plot(mask);
    plot(y_trim);
    legend({'y-trans', 'y-soft thres', 'mask', 'y-trim'});
    % plot(x');
    % plot(y + mean(x));
    % plot(y_trim + mean(x));
end



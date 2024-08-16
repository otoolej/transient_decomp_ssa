%-------------------------------------------------------------------------------
% shorttime_iter_SSA_decomp: Iterative SSA-decomposition using short-time window approach [1]
%
% Syntax: [x_st] = shorttime_iter_SSA_decomp(x, params)
%
% Inputs: 
%     x      - input NIRS signal (1D vector, length-N)
%     fs     - sampling frequency (scalar)
%     params - parameters object 
%
% Outputs: 
%     x_st - structure containing the following fields:
%              x_st.nirs      : input signal 'x''
%              x_st.component : extracted component (1D vector, length-N)
%              x_st.time      : time in seconds (1D vector, length-N)
%
% Example:
%     %  load the default parameters and a test signal:
%     params = decomp_PARAMS;
%     d = load([params.DATA_DIR 'test_signal.mat']);
%     
%     % extend test signal with copy of itself:
%     x2 = [d.x_test d.x_test];
%     fs = 1 / 6;
%     db_plot = true;
%     
%     % short-time, iterative approach to decomposition
%     x_st = shorttime_iter_SSA_decomp(x2, fs, params, db_plot);
%     
%     
%
% [1] O'Toole JM. Dempsey EM, Boylan GB (2018) 'Extracting transients from cerebral
% oxygenation signals of preterm infants: a new singular-spectrum analysis method' in
% Int Conf IEEE Eng Med Biol Society (EMBC), IEEE, pp. 5882--5885 
% https://doi.org/10.1109/EMBC.2018.8513523


% John M. O' Toole, University College Cork
% Started: 09-02-2021
%
% last update: Time-stamp: <2021-04-27 16:51:10 (otoolej)>
%-------------------------------------------------------------------------------
function x_st = shorttime_iter_SSA_decomp(x, fs, params, db_plot)
if(nargin < 3 || isempty(params)), params = decomp_PARAMS; end
if(nargin < 4 || isempty(db_plot)), db_plot = false; end


db_verbose = false;


x = x(:).';
N = length(x);

if(any(isnan(x)))
    error('Remove NaNs first (interpolate?)');
end



% zero mean (or median) signal:
x_mean = nanmedian(x);
x = x - x_mean;



%---------------------------------------------------------------------
% 1. get information on windowing parameters
%    (e.g. may want to change this for different type of window)
%---------------------------------------------------------------------
% use rectangle window with overlap:
[L_hop, L_epoch, win_epoch]=gen_epoch_window(params.st_overlap, ...
                                             params.st_window_length, 'rect', fs);



%---------------------------------------------------------------------
% if does not fit into epochs then extrapolat the end to fit
%---------------------------------------------------------------------
N_epochs = ceil( (N - L_epoch)/L_hop );

N_epochs = N_epochs + 1;
nw = 0:L_epoch - 1;


N_pad = (N_epochs - 1) * L_hop + L_epoch;
if(N_pad <= N)
    N_pad = N; 
else
    % reflection extrapolation:
    x = [x fliplr(x(N - (N_pad - N):end - 1))];
end
ttime = (0:(length(x) - 1)) ./ (1/6);

if(db_verbose)
    fprintf('N_epochs=%d; Window length=%d; hop=%d; N=%d; N_pad=%d\n', ...
            L_epoch, length(win_epoch), L_hop, N, N_pad);    
end



z_all = zeros(1, N_pad); 
win_summed = zeros(1, N_pad);
x_component = zeros(1, N_pad); 


for k = 0:N_epochs - 1
    nf = mod(nw + k * L_hop, N_pad);
    x_epoch = x(nf + 1) .* win_epoch.';

    DBplot_iter = false;
    if(DBplot_iter)
        set_figure(99);
        subplot(312); 
        plot(ttime(nf + 1), x_epoch);
        subplot(311); 
        plot(ttime, x);
        ys = ylim;
        hrec = rectangle('position',[ttime(nf(1) + 1) ys(1) ttime(nf(end) + 1) - ...
                            ttime(nf(1) + 1) ys(2) - ys(1)], ...
                         'facecolor',[1 1 1] .* 0.8);
        uistack(hrec,'bottom');
    end    
    
    
    %---------------------------------------------------------------------
    % SSA
    %---------------------------------------------------------------------
    x_ssa_dct = iterative_SSA_decomposition(x_epoch, params.L_ssa_ev, params.SSA_METHOD, ...
                                            params.ITER_L_ssa_ev, params.USE_DCT, ...
                                            params.DCT_CUTOFF, false);
    
    
    %---------------------------------------------------------------------
    % post processing
    %---------------------------------------------------------------------
    x_ssa_dct = trim_transient_component(x_ssa_dct, fs, false);

    
    if(DBplot_iter)
        set_figure(99, 0);
        subplot(313); plot(ttime(nf + 1), x_ssa_dct);
        disp('--- paused; hit key to continue ---'); pause;
    end

    
    x_component(nf + 1) = x_component(nf + 1) + x_ssa_dct;
    z_all(nf + 1) = z_all(nf + 1) + x_epoch; 
    win_summed(nf + 1) = win_summed(nf + 1) + win_epoch.';            
end


%---------------------------------------------------------------------
% reconstruct component and residual signal
%---------------------------------------------------------------------
z_all = z_all ./ win_summed;
x_component = x_component ./ win_summed;
% add back mean:
z_all = z_all + x_mean;
x_component = x_component + x_mean;
x = x + x_mean;

% remove the extrapolated parts:
z_all = z_all(1:N); 
x = x(1:N);
x_component = x_component(1:N);
ttime = ttime(1:N);


% put in structure:
x_st.nirs = x;
x_st.component = x_component;
x_st.time = ttime;


%---------------------------------------------------------------------
% plot?
%---------------------------------------------------------------------
if(db_plot)
    set_figure(1); 
    ttime = x_st.time / 3600;
    
    hx(1) = subplot(3, 1, 1);
    plot(ttime, x_st.nirs);
    title('NIRS signal');
    hx(2) = subplot(3, 1, 2);
    plot(ttime, x_st.component);
    title('extracted component');
    hx(3) = subplot(3, 1, 3);
    hold all;
    plot(ttime, x_st.nirs);
    plot(ttime, x_st.nirs - x_st.component + x_mean);
    title('signal minus component');
    legend({'signal', 'signal - component'});
    linkaxes(hx, 'x');
end


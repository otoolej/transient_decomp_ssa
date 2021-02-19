%-------------------------------------------------------------------------------
% iterative_SSA_decomposition: iteratative process of SSA decomposition from [1]
%
% Syntax: y_ssa=iterative_SSA_decomposition(x,L,method,all_Ls,use_dct,dct_keep,db_plot)
%
% Inputs: 
%     x        - input signal (1D vector, length-N)
%     L        - embedding dimension (labelled 'M' in [1])
%     method   - method to determine threshold between signal/noise components
%                either: 'vautard_ghil', 'celka_mdl', or 'eigenvalue_ranking'
%     all_Ls   - array of embedding dimensions
%     use_dct  - use the discrete cosine transform (true/false)
%     dct_keep - proportion of DCT transform to keep
%     db_plot  - plot (true/false) 
%
% Outputs: 
%     
%     % load the default parameters and a test signal:
%     params = decomp_PARAMS;
%     d = load([params.DATA_DIR 'test_signal.mat']);
%     
%     % set parameters (default values):
%     L = params.L_ssa_ev;
%     method = params.SSA_METHOD;
%     iter_Ls = params.ITER_L_ssa_ev;
%     use_dct = params.USE_DCT;
%     frac_dct_ignore = params.DCT_CUTOFF;
%     db_plot = true;
%     
%     % extract the transients:
%     y_ssa = iterative_SSA_decomposition(d.x_test, L, method, iter_Ls, use_dct, frac_dct_ignore, db_plot);
%     
%
% Example:
%     
%
% [1] O'Toole JM. Dempsey EM, Boylan GB (2018) 'Extracting transients from cerebral
% oxygenation signals of preterm infants: a new singular-spectrum analysis method' in
% Int Conf IEEE Eng Med Biol Society (EMBC), IEEE, pp. 5882--5885 
% https://doi.org/10.1109/EMBC.2018.8513523

% John M. O' Toole, University College Cork
% Started: 11-12-2017
%
% last update: Time-stamp: <2021-02-16 15:28:48 (otoolej)>
%-------------------------------------------------------------------------------
function [y_ssa]=iterative_SSA_decomposition(x,L,method,all_Ls,use_dct,dct_keep,db_plot)
if(nargin<2 || isempty(L)), L=20; end
if(nargin<3 || isempty(method)), method='Vautard_Ghil'; end
if(nargin<4 || isempty(all_Ls)), all_Ls=[25:5:40]; end
if(nargin<4 || isempty(use_dct)), use_dct=1; end
if(nargin<5 || isempty(dct_keep)), dct_keep=0.1; end
if(nargin<7 || isempty(db_plot)), db_plot=0; end

DBverbose=0;

% important to keep mean of original signal:
xmean=nanmean(x);

%---------------------------------------------------------------------
% first pass:
%---------------------------------------------------------------------
subtract_mean=true;
[y_decomp,S]=noise_reduction_SSA(x,L,method,use_dct,dct_keep,subtract_mean,[],db_plot); 

% if no iterative component:
if(isempty(all_Ls))
    y_ssa=y_decomp;
    return;
end


%---------------------------------------------------------------------
% iterative over all L
%---------------------------------------------------------------------
y_iter=y_decomp; 
frac_LS=S/L;
subtract_mean=false;
for n=1:length(all_Ls)
    
    % set limit of max. S
    L_max=round(all_Ls(n)*frac_LS);

    if(DBverbose)
        fprintf('L=%d (L-max=%d)\n',all_Ls(n),L_max);
    end
    
    
    ymean=nanmean(y_iter);

    [y_iter,S]=noise_reduction_SSA(y_iter-ymean,all_Ls(n),method, ...
                                   use_dct,1,subtract_mean,L_max,false); 
    
    
    frac_LS=S/all_Ls(n);
    
    % important, sum back original mean:
    y_iter=y_iter+xmean;
end


y_ssa=y_iter;


db_plot = false;
if(db_plot)
    set_figure(44);
    subplot(2,1,1); hold all;
    plot(x); plot(y_ssa); 
    legend({'x + transient', 'iter. SSA decomp.'});
    subplot(2,1,2); hold all;
    plot(x); plot(x - y_ssa + xmean);
    legend({'x', 'x - iter. SSA decomp.'});
end

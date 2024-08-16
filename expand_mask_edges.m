%-------------------------------------------------------------------------------
% expand_mask_edges: expand the edges of binary mask
%
% Syntax: y = expand_mask_edges(x)
%
% Inputs: 
%     x - discrete mask E [-1, 0, 1]
%
% Outputs: 
%     y - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 22-04-2021
%
% last update: Time-stamp: <2021-04-22 17:18:41 (otoolej)>
%-------------------------------------------------------------------------------
function y = expand_mask_edges(x, L_win)

N = length(x);
w = blackmanharris(L_win);
wh = w(1:floor(length(w) / 2));
Lh = length(wh);
    
% a) first the +ve spikes:
yp = x;
yp(yp < 0) = 0;
yp = smooth_edges(yp, wh, Lh, N);

% b) then the -ve spikes:
yn = x;
yn(yn > 0) = 0;
yn = smooth_edges(-yn, wh, Lh, N);

% combine and hard-limit:
y = yp + yn;
y(y > 1) = 1;



function x = smooth_edges(x, w, L, N)
%---------------------------------------------------------------------
% smooth edges with the window
%---------------------------------------------------------------------
[lens, istart, iend] = len_cont_zeros(x, 1);

for n = 1:length(istart)
    idxs = (istart(n) - L):(istart(n) - 1);
    idxe = (iend(n) + 1):(iend(n) + L);    
    
    % make sure are within limits:
    idxs(idxs < 1) = 1;
    idxs(idxs > N) = N;    
    idxe(idxe < 1) = 1;
    idxe(idxe > N) = N;    

    x(idxs) = x(idxs) + w';
    x(idxe) = x(idxe) + w(end:-1:1)';
end

x(x > 1) = 1;

    

    
    

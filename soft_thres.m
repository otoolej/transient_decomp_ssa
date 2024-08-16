%-------------------------------------------------------------------------------
% soft_thres: soft thresholding
%
% Syntax: y = soft_thres(x, L)
%
% Inputs: 
%     x, L - 
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
% last update: Time-stamp: <2021-04-22 19:10:26 (otoolej)>
%-------------------------------------------------------------------------------
function y = soft_thres(x, L)
    
y = zeros(size(x));
ithres = find(abs(x) > L);
y(ithres) = x(ithres) - L * sign(x(ithres));


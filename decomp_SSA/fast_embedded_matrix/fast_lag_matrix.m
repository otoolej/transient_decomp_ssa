%#codegen
%-------------------------------------------------------------------------------
% fast_lag_matrix: implement lagged matrix
%
% Syntax: x_lag=fast_lag_matrix(x,L)
%
% Inputs: 
%     x,L - 
%
% Outputs: 
%     x_lag - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 05-12-2017
%
% last update: Time-stamp: <2017-12-05 11:19:09 (otoolej)>
%-------------------------------------------------------------------------------
function x_lag=fast_lag_matrix(x,L)

N=int32(length(x));

M=N-L+1;
L_1=L-1;


x_lag=zeros(L,M);


%---------------------------------------------------------------------
% 1. lag matrix
%---------------------------------------------------------------------
for k=1:M
    x_lag(:,k)=x( k:(k+L_1) );
end



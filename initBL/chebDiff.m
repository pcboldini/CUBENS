% Routine for the calculation of Chebychev differentiation matrices from:
% Lloyd N. Trefethen, Spectral Methods in MATLAB, SIAM, Philadelphia, 2000
% http://people.maths.ox.ac.uk/trefethen/cheb.m
%
%

function [D] = chebDiff(N)

if N==0, D=0; return, end
N=N-1;
x = -cos(pi*(0:N)/N)'; 
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';                  
D  = (c*(1./c)')./(dX+(eye(N+1)));
D  = D - diag(sum(D'));


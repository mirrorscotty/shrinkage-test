% INVLAP  Numerical Inversion of Laplace Transforms 
function [s]=INVLAPs(tini,tend,nnt,a,ns,nd);
% tini, tend are limits of the solution interval
% nnt is total number of time instants
% a, ns, nd are parameters of the method 
% if not given, the method uses implicit values a=6, ns=20, nd=19
% it is recommended to preserve a=6
% increasing ns and nd leads to lower error
% an example of function calling  
if nargin==3
  a=6; ns=20; nd=19;  end;    % implicit parameters
radt=linspace(tini,tend,nnt); % time vector
if tini==0  radt=radt(2:1:nnt);  end;  % t=0 is not allowed
%tic					% measure the CPU time
for n=1:ns+1+nd               % prepare necessary coefficients
   alfa(n)=a+(n-1)*pi*j;
end;
for kt=1:nnt                  % cycle for time t
   tt=radt(kt);
   s=alfa/tt;                 % complex frequency s
end;
%toc
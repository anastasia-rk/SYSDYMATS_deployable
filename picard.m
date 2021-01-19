function eta = picard(U,s,b,d)
%PICARD Visual inspection of the Picard condition - from Mathworks
%
% eta = picard(U,s,b)
% eta = picard(U,s,b,d)
% eta = picard(U,sm,b)    ,  sm = [sigma,mu]
% eta = picard(U,sm,b,d)
%
%  eta(i) = |U(:,i)'*b|/s(i).  U and s are from compact svd
%
% The smoothing is a geometric mean over 2*d+1 points centered at point #i.
% If nargin = 3, then d = 0 (i.e, no smothing).

% Reference: P. C. Hansen, "The discrete Picard condition for discrete
% ill-posed problems", BIT 30 (1990), 658-672.


% Initialize
[n,ps] = size(s); beta = abs(U(:,1:n)'*b); eta = zeros(n,1);
if (nargin==3), d = 0; end
if (ps==2), s = s(:,1)./s(:,2); end
d21 = 2*d+1; nEta = 1+d:n-d;
if ~all(s), warning('Division by zero singular values'), end
w = warning('off');


for i=nEta
  eta(i) = (prod(beta(i-d:i+d))^(1/d21))/s(i);
end
warning(w);

% Plot the data.
semilogy(1:n,s,'o',1:n,beta,'x',nEta,eta(nEta),'.-')
xlabel('i');ylabel('SVD decompostion')
if (ps==1)
  legend('$\sigma_i$','$|u_i^T Y|$','$|u_i^T Y|/\sigma_i$')
else
  legend('$\sigma_i/\mu_i$','$|u_i^Tb|$','$|u_i^Tb| / (\sigma_i/\mu_i)$',...
      'Location','NorthWest')
end
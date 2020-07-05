function [r_km] = residual_update(r_km_1,q_km)
% Notation mataches that of (Wei,Lang,Billings, 2008) and differs from the main file
% Input
%   r_km_1  - residual  for the previous regression model;
%   q_km    - orthogonal basis vector;
% Output
%   r_km - the new residual
    coeff = (r_km_1'*q_km)/(q_km'*q_km);
    r_km  = r_km_1 - coeff*q_km;
end
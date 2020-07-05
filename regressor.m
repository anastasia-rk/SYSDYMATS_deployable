function [result] = regressor(x,indeces)
% Input
%   x - input vector x = [y(t-n_y:t-1); u(t-n_u:t-1)];
%   indeces - elements of x present in the polynomial term;
% Output
%   result - the product
x_ind = x(indeces,:); % select the 
result = prod(x_ind);
end
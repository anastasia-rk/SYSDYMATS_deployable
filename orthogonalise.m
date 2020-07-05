function [p_new] = orthogonalise(phi_kj,q_prev,iTerm)
% Notation mataches that of (Wei,Lang,Billings, 2008) and differs from the main file
% Input
%   phi_kj  - basis vector (size nNarx x 1);
%   q_prev - orthogonalised basis vectors from iteration 0 to iTerm-1 (size (nNarx x iTerm));
% Output
%   p_new - the new basis corresponding to orthogonalised significant vectors 
    sum_q = zeros(size(phi_kj));
    for s=1:iTerm-1                                                         % Over all significant terms
        q = q_prev(:,s);                                                    % Get the orthogonalised basis vector
        sum_q = sum_q + phi_kj'*q/(q'*q)*q;                                 % sum all orthogonal terms
    end
    p_new  = phi_kj - sum_q;                                                % new basis vector 
end
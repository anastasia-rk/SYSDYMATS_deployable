function [Q, R] = rdqr(A,tol)
    % Find the QR factorization with pivoting
    [Q, R] = mgson(A);

    % Since the diagonal of Q has elements decreasing in magnitude, find the
    % first element that's near zero
    rrank = find(abs(diag(R))<tol);
    % If we have zero rank, throw an error
    if isempty(rrank)
        reduced = rank(R);
    else
    if rrank<=1
        error('Attempted to take a QR factorization of a rank-0 matrix')
    end
        reduced = rrank(1)-1;
    end

    % Reduce the size of Q and R accordingly
    Q = Q(:,1:reduced);
    R = R(1:reduced,:);
end
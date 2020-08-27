% Create a tall and skinny, rank deficient matrix
m = 10;
n = 5;
r = 3;
A = randn(m,r)*randn(r,n);

% Find the rank deficient QR factorization with pivoting
[Q R p] = qrrd(A,1e-12);

% Find the residual
fprintf('(tall) || A(:,p) - QR ||: %e\n',norm(A(:,p)-Q*R,'fro'));

% Make A short and fat by taking the transpose
A = A';

% Repeat the exercise
[Q R p] = qrrd(A,1e-12);
fprintf('(short) || A(:,p) - QR ||: %e\n',norm(A(:,p)-Q*R,'fro'));


% Finds a QR factorization of a rank deficient matrix
function [Q R p] = qrrd(A,tol)
    % Find the QR factorization with pivoting
    [Q R p] = qr(A,0);

    % Since the diagonal of Q has elements decreasing in magnitude, find the
    % first element that's near zero
    rrank = find(abs(diag(R))<tol);
    rrank = rrank(1)-1;

    % If we have zero rank, throw an error
    if rrank==0
        error('Attempted to take a QR factorization of a rank-0 matrix')
    end

    % Reduce the size of Q and R accordingly
    Q = Q(:,1:rrank);
    R = R(1:rrank,:);
end
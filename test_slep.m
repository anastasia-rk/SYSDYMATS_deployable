% clear, clc;
Y = Y_train{1}(1:35000,:);
A1 = M_train{1}(1:35000,:);
A_tilda = [];
ind_orig = [];
for i = 1:max(p_sparesgroup)
    index = find(p_sparesgroup == i);
   A_tilda = [A_tilda A1(:,index)];
   ind_orig = [ind_orig; index];
end

%----------------------- Set optional items -----------------------
opts=[];

% Starting point
opts.init=2;        % starting from a zero point

% Termination 
opts.tFlag=5;       % run .maxIter iterations
opts.maxIter=100;   % maximum number of iterations

% regularization
opts.rFlag=1;       % use ratio

% Normalization
opts.nFlag=0;       % without normalization

% Group Property
opts.ind=[ [1, 11, sqrt(11)]', [12, 22, sqrt(11)]',...
    [23, 33, sqrt(11)]', [34, 44, sqrt(11)]', [45,55, sqrt(11)]',  [56,66, sqrt(11)]'];

%----------------------- Run the code -----------------------
alpha = 0.5;
lambda = 0.001;
z=[lambda*(1-alpha),lambda*alpha];
tic;
[x1, funVal1, ValueL]= sgLeastR(A_tilda, Y, z, opts);
toc;

x_new = x1(ind_orig);
Betas_spgl = reshape(x_new,[finalTerm L])
% b=x1~=0;
% nnz(b)
% 
y_new = A_tilda*x1;
figure;
plot(Y(1:200)); hold on
plot(y_new(1:200))

%%regroup terms again accoring to the original model;


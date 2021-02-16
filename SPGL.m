function [X] = SPGL(Y,A,lambda,alpha,groupIndex,X_lasso)
% modify regressor matrix to order the groups
A_tilda     = [];
X_0         = [];
ind_orig    = [];
group_info  = [];
nGroups     = max(groupIndex);
ind_end     = 0;
for i = 1:nGroups
   index = find(groupIndex == i);
   ind_start = ind_end + 1; % group beginning index
   ind_end = ind_end + length(index); % group ending index
   A_tilda = [A_tilda A(:,index)]; % re-order regressors
   X_0     = [X_0; X_lasso(index)];
   ind_orig = [ind_orig; index]; % memorise original positions
   new_group = [ind_start, ind_end, sqrt(length(index))]'; % form new group info
   group_info = [group_info new_group]; 
end
L = length(index);
% beta_0 = M_tilda\Y;
%----------------------- Set optional items -----------------------
opts=[];

% Starting point
opts.init = 1;        % starting from the LS estimate
opts.x0 = X_0;
% Termination 
opts.tFlag=4;       % run until norm of residual below threshold
% opts.maxIter=200;   % maximum number of iterations

% regularization
opts.rFlag=1;       % use ratio

% Normalization
opts.nFlag=0;       % without normalization

% Group Property - define according to the group structure
opts.ind=group_info;

%----------------------- Run the code -----------------------
z=[lambda*(1-alpha),lambda*(alpha)];
[x1, funVal1, ValueL]= sgLeastR(A_tilda, Y, z, opts);
X = reshape(x1,[L,nGroups])';
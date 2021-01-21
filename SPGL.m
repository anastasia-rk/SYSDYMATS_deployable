function [beta] = SPGL(Y,M,lambda,alpha,groupIndex)
% modify regressor matrix to order the groups
M_tilda = [];
ind_orig = [];
group_info = [];
nGroups = max(groupIndex);
ind_end   = 0;
for i = 1:nGroups
   index = find(groupIndex == i);
   ind_start = ind_end + 1; % group beginning index
   ind_end = ind_end + length(index); % group ending index
   M_tilda = [M_tilda M(:,index)]; % re-order regressors
   ind_orig = [ind_orig; index]; % memorise original positions
   new_group = [ind_start, ind_end, sqrt(length(index))]'; % form new group info
   group_info = [group_info new_group]; 
end
beta_0 = M_tilda\Y;
%----------------------- Set optional items -----------------------
opts=[];

% Starting point
opts.init = 2;        % starting from the LS estimate
% opts.x0 = beta_0;
% Termination 
opts.tFlag=4;       % run until norm of residual below threshold
% opts.maxIter=100;   % maximum number of iterations

% regularization
opts.rFlag=0;       % use ratio

% Normalization
opts.nFlag=0;       % without normalization

% Group Property - define according to the group structure
opts.ind=group_info;

%----------------------- Run the code -----------------------
z=[lambda*(1-alpha),lambda*alpha];
[b1, funVal1, ValueL]= sgLeastR(M_tilda, Y, z, opts);
beta = b1(ind_orig);

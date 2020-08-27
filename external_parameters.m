local_init;
%% dataset C,D
% Dictionary
params_all = {'D_init','D_imposed','L_init','L_imposed','R_imposed',...
         'A_imposed','V_imposed','D_relax','L_cut','W_relax','Dens_relax'};
L = numel(params_all);
values(:,1)  = [30,30,30,30,30,50,50,50,50,50]';
values(:,2)  = [19,19,19,19,19,19,19,19,19,19]';
values(:,3)  = [170,150,130,110,90,160,140,120,100,90]';
values(:,4)  = [80,80,80,80,80,75,75,75,75,75]';
values(:,5)  = [1.58,1.58,1.58,1.58,1.58,2.63,2.63,2.63,2.63,2.63]';
values(:,6)  = [2.13,1.88,1.63,1.38,1.13,2.13,1.87,1.60,1.33,1.20]';
values(:,7)  = [5.30,4.67,4.05,3.43,2.80,14.8,12.9,11.1,9.2,8.3]';
values(:,8)  = [20,20,20,20,20,20,20,20,20,20]';
values(:,9)  = [74,72,69,61,56,67,66,64,62,57]';
values(:,10) = [2.87,2.56,2.26,1.83,1.64,5.8,5.28,4.48,4.2,3.45]';
values(:,11) = [0.123,0.113,0.104,0.095,0.093,0.276,0.255,0.230,0.209,0.193]';
%% total number of permutations
M       = permn(1:L,2);                                                     % get all permutations with repetition
ind     = find(M(:,2)>=M(:,1));                                             % sort out only increasing indeces
indeces = M(ind,:);                                                         % updated set
K       = length(indeces)                                                   % total number of sample pairs
clear ind
%% Check correlations of sample pairs
check_var = 0/0;
% for all pairs
for k = 1:K   
    % perform correlation test (possible types: Kendall, Soearman, Pearson)
    [rho(k),pval(k)]= corr(values(:,indeces(k,1)),values(:,indeces(k,2)));
end
% Find all pairs that don't reject non-correlation hypothesis
ind          = find(~isnan(pval) & pval>0.0005);
% update arrays
rho          = rho(ind);
pval         = pval(ind);    
indeces_new  = indeces(ind,:);
[rho_max,k_max] = min(rho);
param_select = params_all(indeces_new(k_max,:));

my_choice = [9 11];
values = values(:, my_choice);
params{1} = params_all{1,my_choice(1)};
params{2} = params_all{1,my_choice(2)};
fileName = 'External_parameters_C';
save(fileName,'params','values');
fileName = 'External_parameters_D';
save(fileName,'params','values');
clear values params
%% Dataset S-Z
params = {'Cut length, m'};
values = [3 3 3 5 5 5 7 7 7 0 0 0]';
values = values./1000;
% values = [values_cut; values_cut; values_cut; values_cut]
fileName = 'External_parameters_S';
save(fileName,'params','values');
fileName = 'External_parameters_Z';
save(fileName,'params','values');
fileName = 'External_parameters_Y';
save(fileName,'params','values');

% Tests below have been done five times for each specimen with varying
% noise power
I = ones(1,5);
%% Datasets HF - hardfoam
params = {'Thickness, m','Tension'};
values1 = [5*I 5*I 10*I 10*I 10*I 10*I 10*I 15*I 15*I 5*I 10*I 15*I]./1000;...
values2 = [0*I 20*I 0*I 5*I 10*I 15*I 20*I 0*I 20*I 0*I 0*I 0*I]./100;
values  = [values1; values2]';
fileName = 'External_parameters_HF';
save(fileName,'params','values');

%% Daatsets H,V - soft foam
params = {'Compression ratio'};
values = [40*I 50*I 60*I 70*I 80*I 0*I]';
values = values./100;
fileName = 'External_parameters_VS';
save(fileName,'params','values');
fileName = 'External_parameters_HS';
save(fileName,'params','values');






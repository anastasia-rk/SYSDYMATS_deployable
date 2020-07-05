function main(foamset,dataset,regressors,direction)
local_init;
%%
switch foamset
    case 'foam_2010'
            Files = [1 2 4 5 6 7 9 10];
            testFiles  = [3 8];
    case 'foam_2019'
            Files = [3 6 9 12];
            testFiles  = [1 3 5 7 8 10 11];
    case 'foam_2020'
    switch dataset
        case 'HS'
            Files = [4 5 8 10 13 15 18 20 22 25];
            testFiles  = [1 3 8 9 12 14 17 23];
        case 'VS'
            Files = [4 5 8 10 13 15 18 20 22 25]; 
            testFiles  = [1 3 8 9 12 14 17 23];
        case 'HF'
            Files = [1 10 15 20 30 35 40 41]; 
            testFiles  = [5 8 11 19 25 29 37 45];
    end 
    case 'sweep_2020'
         Files = [1 2 3 5];
         testFiles  = [1 3 4 6];       
end
folder = 'results_cv';                                                     % specify category where to save files
dFolder = 'dictionaries';
switch regressors
    case 'Shift'
        metaFileName = ['Meta_',dataset];
        load(metaFileName);
        names = {'set','ny','nu'};                                          % names used to define results folder name (no more than 3).
        if normC ~= 1
            folder = [folder,'_norm'];
            dFolder = [dFolder,'_norm'];
        end
        folderName = make_folder(folder,names,dataset,n_y,n_u);             % create results folder
        dictFolder = make_folder(dFolder,names,dataset,n_y,n_u);            % create results folder
        d       = n_y + n_u;                                                % size of input vector x
    case 'Delta'
         folder = ['delta_',folder];
         dFolder = ['delta_',dFolder];
        switch direction
            case 'Backward'
                folder = [folder,'_b'];
                dFolder = [dFolder,'_b'];
                 metaFileName = ['Meta_delta_',dataset,'_b'];
            case 'Forward'
                folder = [folder,'_f'];
                dFolder = [dFolder,'_f'];
                 metaFileName = ['Meta_delta_',dataset,'_f'];
        end
        load(metaFileName);
        if normC ~= 1
            folder = [folder,'_norm'];
            dFolder = [dFolder,'_norm'];
        end
         names = {'set','lambda'};
         folderName = make_folder(folder,names,dataset,lambda);             % create results folder
         dictFolder = make_folder(dFolder,names,dataset,lambda);            % create results folder
         n_u = 1+lambda;
         n_y = 1;
         d = lambda*2;
end
dict_set = ['dict_',dataset];                                   
fileNames = sym(dict_set,[1 K]);                                            % vector of filenames
K = length(Files);
% Set maximum number of covariates
if length(dict_terms) < 30                                                  % Maximum significant terms (if the algorithm is not terminated by the criterion)
    maxSign = length(dict_terms);
else
    maxSign = 30;                                                           
end
dict_terms_all = dict_terms;
%% Create the regression matrix based on the dataset (does not depend on CV parameters)
load(['External_parameters_',dataset]);
x = values(Files,1);                                                        % design parameter values for training
x_valid = values(testFiles,1);                                              % design parameter values for validation
if size(values,2) > 1
    y = values(Files,2);
    y_valid = values(testFiles,2);
else 
    y = [];
end
A = ones(size(x));                                                          % create unit vector for constants in design matrix
A_valid = ones(size(x_valid));                                              % create unit vector for constants in validation matrix 
A_symb{1} = '1';
if ~isempty(y)                                                              % unknown mapping is a surface
   powers = permn(0:2,2);                                                   % permuntations of all 
   powers = powers(2:end,:);    
   nCols = min(size(powers,1),K);                                           % number of terms in the model shouldn't be higher then K
   for iCol = 1:nCols
       xCol = x.^powers(iCol,1);
       yCol = y.^powers(iCol,2);
       A = [A xCol.*yCol];
       xvCol = x_valid.^powers(iCol,1);
       yvCol = y_valid.^powers(iCol,2);
       A_valid = [A_valid xvCol.*yvCol];
       A_symb{iCol+1} = ['$x^',num2str(powers(iCol,1)),'$ $y^',num2str(powers(iCol,2)),'$']; 
   end
else                                                                        % unknown mapping is a curve
    nCols =  min(3,K);                                                      % number of terms in the model shouldn't be higher then K
    for iCol = 1:nCols                                                      % limit order of the model by the number of experimants
       A = [A x.^(iCol)];
       A_valid = [A_valid x_valid.^(iCol)];
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Structure identification with CV model selection
    index       = 1:nNarx;                                                  % structure id over all times
    dict_terms  = dict_terms_all;                                           % refill the dictionary
    %% Select first significant basis vector for all datasets
    iTerm = 1;                                                              % the first significant term
    AERR{iTerm} = zeros(nTerms,1);                                          % placeholder for AERR criteria
    for iFile = Files                                                       % over all datasets
        fName   = [dictFolder,'/',char(fileNames(iFile))];
        File    = matfile(fName,'Writable',true);
        y_n     = File.y_narx(:,1);
        residual_init{iFile} =  y_n(index,1);                               % initial residual
        for jTerm = dict_terms                                              % over all polynomial terms in the dictionary
            term_all    = File.term(:,jTerm);
            term0       = term_all(index,:);
            cf(iFile,jTerm)     = cor_sqr(residual_init{iFile},term0);      % squared correlation coefficient for the dataset and the polynomial term
            AERR{iTerm}(jTerm)  = AERR{iTerm}(jTerm) + cf(iFile,jTerm);     % Average error reduction ration over all datasets
            clear term0 term_all
        end
        clear File y_n 
    end
    AERR{iTerm}(:,:)    = AERR{iTerm}(:,:)/K;
    [AERR_m,iMax]       = max(AERR{iTerm});                                 % find the index of the term with the highest criterion across all datasets
    AERR_mm(iTerm,1)    = AERR_m;
    S(iTerm)            = iMax;                                             % save index of the term
    dict_terms(iMax)    = [];                                               % reduce the dictionary of available term
    BIC_sum             = 0;
    for iFile = Files                                                       % over all datasets
        fName = [dictFolder,'/',char(fileNames(iFile))];
        File  = matfile(fName,'Writable',true);
        term_all    = File.term(:,iMax);
        alpha{iFile}(:,iTerm)    = term_all(index,:);                       % the corresponding basis candidate term    
        phi  {iFile}(:,iTerm)    = term_all(index,:);                       % the corresponding basis vector 
        residual{iFile}(:,iTerm) = residual_update(residual_init{iFile},... % the corresponding model residual
                                               phi{iFile}(:,iTerm));                                                        
        BIC_sum  = BIC_sum  +  BIC(residual{iFile}(:,iTerm),nNarx,iTerm);   % BIC for the iFile dataset
        clear File term_all
    end
    BIC_all(iTerm)            = BIC_sum/K;                                  % average AMDL over all sets
    significant_term{iTerm}   = symb_term{S(iTerm)};
%% Main loop   
    converged   = false;
    bics        = [];
    while(iTerm < maxSign) %&& ~converged                                   % loop over the number of significant terms
        iTerm = iTerm + 1;                                                  % increase the number of significant terms
        AERR{iTerm} = zeros(nTerms,1);                                      % placeholder for AERR criteria
        for iFile = Files                                                   % over all datasets
            fName   = [dictFolder,'/',char(fileNames(iFile))];
            File    = matfile(fName,'Writable',true);
            for jTerm = dict_terms                                          % over all polynomial terms in the dictionary
                term_all    = File.term(:,jTerm);
                p{iTerm,iFile}(:,jTerm) = orthogonalise(term_all(index,:),...
                                                    phi{iFile},iTerm);      % orthogonalise basis
                cf(iFile,jTerm)         = cor_sqr(residual_init{iFile},...
                                              p{iTerm,iFile}(:,jTerm));     % squared correlation coefficient for the dataset and the polynomial term
                AERR{iTerm}(jTerm) = AERR{iTerm}(jTerm) + cf(iFile,jTerm);  % average error reduction ration over all datasets
            end
            clear File
        end
        AERR{iTerm}(:,:)    = AERR{iTerm}(:,:)/K;
        [AERR_m,iMax]       = max(AERR{iTerm});                             % find the index of the term with the highest criterion across all datasets
        AERR_mm(iTerm,1)    = AERR_m;
        S(iTerm)            = iMax;                                         % save index of the term  
        ind = find(dict_terms == iMax);
        dict_terms(ind) = [];                                               % Reduce the dictionary of available terms
        BIC_sum         = 0;
        for iFile = Files
            fName   = [dictFolder,'/',char(fileNames(iFile))];
            File    = matfile(fName,'Writable',true);
            alpha{iFile}(:,iTerm) = File.term(index,S(iTerm));              % the corresponding basis candidate term    
            phi{iFile}(:,iTerm)   = p{iTerm,iFile}(:,S(iTerm));             % the corresponding basis vector 
            residual{iFile}(:,iTerm) = residual_update(residual{iFile}(:,iTerm-1),...
                                                   phi{iFile}(:,iTerm));    % the corresponding model residual                                 
            BIC_sum  = BIC_sum  +  BIC(residual{iFile}(:,iTerm),nNarx,iTerm); % BIC for the iFile dataset
            clear File x_n
        end
        significant_term{iTerm} = symb_term{S(iTerm)};
        BIC_all(iTerm) = BIC_sum/K;                                         % average AMDL over all sets
        converged_BIC = (abs((BIC_all(iTerm) - BIC_all(iTerm-1))/BIC_all(iTerm)) < 0.0002); % check convergence
        if converged_BIC
            bics  = [bics,iTerm];
        end
    end
if isempty(bics)
    finalTerm = 20;
else
    finalTerm = bics(1);
end
BIC_trunc = BIC_all(1:finalTerm);
fig('BIC',visFlag);
plot([1:maxSign],BIC_all); hold on;
plot(finalTerm,BIC_trunc(end),'*');
xlabel('Terms');ylabel('BIC');
tikzName = [folderName,'/BIC_all_folds.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '6cm', 'width','6cm','checkForUpdates',false); 
%% Create column of names
    for iTerm=1:finalTerm
        temp = arrayfun(@char, significant_term{iTerm}, 'uniform', 0);
        f_mpo{iTerm} = inline(inline_term{S(iTerm)},'u','y','t');
        if length(temp) > 0
            str = temp{1};
            for iString=2:length(temp)
                str = [str,temp{iString}];
            end
        end
        Terms{iTerm,1} = strcat('$',str,'$');
        clear temp
    end
    Step = [1:finalTerm]';
    Tab = table(Step,Terms);
    clear AERR
    AERR  = round(AERR_mm(1:finalTerm,1)*100,3);
%% Parameter estimation
    for iFile=Files
        U{iFile} = zeros(finalTerm,finalTerm);                              % placeholder for upper-trig unit matrix
        iTerm = 1;                                                          % for the first term
        g{iFile}(iTerm) = (residual_init{iFile}'*phi{iFile}(:,iTerm))/...
                          (phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));       % top raw of the rotation matrix
        for jTerm =iTerm:finalTerm
            U{iFile}(iTerm,jTerm) = alpha{iFile}(:,jTerm)'*phi{iFile}(:,iTerm)/...
                                (phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));
        end
        for iTerm = 2:finalTerm                                             % loop over significant terms
            g{iFile}(iTerm,1) = (residual{iFile}(:,iTerm-1)'*phi{iFile}...
                    (:,iTerm))/(phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));  % righthand side (normalised)
            for jTerm =iTerm:finalTerm
                U{iFile}(iTerm,jTerm) = alpha{iFile}(:,jTerm)'*phi{iFile}...
                     (:,iTerm)/(phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));  % upper triangular unit matrix
            end
        end
        Dets = det(U{iFile});
        Theta(:,iFile) = linsolve(U{iFile},g{iFile},struct('UT', true));    % solve upper triangular system via backward substitution
        Parameters = round(Theta(:,iFile),2);
        varName = [dataset,num2str(iFile)];
        Tab = addvars(Tab,Parameters,'NewVariableNames',varName);
    end
%% Save internal parameters to table
Tab   = addvars(Tab,AERR,'NewVariableNames',{'AERR($\%$)'});
internalParams = addvars(Tab,BIC_trunc','NewVariableNames',{'BIC'});
if disFlag 
    internalParams
end
tableName = [folderName,'/Thetas_overall'];
table2latex(internalParams,tableName);
clear Tab tableName
clear AERR alpha U phi residual p Theta g
%% Form regression matrices for all time points
% vectorisation for joint estimation
I   = eye(finalTerm);                                                       % unit matrix, size NxN
Kr  = kron(A,I);
L   = size(A,2);
p_sparesgroup = [];
for iBeta = 1:L                                                             % across all design parameter values
    p_sparesgroup = [p_sparesgroup 1:finalTerm];                            % form groups for CVX sparse group lasso
end
p_sparesgroup = p_sparesgroup';
Phi_bar = [];
Y_all   = [];
for iFile = Files
    fName     = [dictFolder,'/',char(fileNames(iFile))];
    File      = matfile(fName,'Writable',true);
    indSign   = S(1:finalTerm);                                             % select the indeces of significant terms from the ordered set
    Phi_file  = File.term(:,:);                                             % extract all terms into a vector - cannot reorder directly in files
    y_file    = File.y_narx(:,:);                                           % extract output
    Phi       = Phi_file(:,indSign);                                        % select only signficant terms
    Phi_bar   = blkdiag(Phi_bar,Phi);                                       % diagonal block, size TKxNK
    Y_all     = [Y_all; y_file];                                            % block vector, size TKx1
end
M_all  = Phi_bar*Kr;  
[Q_all,R_all] = mgson(M_all);                                               % Orthogonalise M via modified Gram-Schmidt
% !!!!!!!!!!!! From this point onwards we work with Q_all for
% cross-validation. We will utilise the R_all for orthogonal solution only after
% optimal regularisation parameters have been found
%% Block CV with K-folds
nFolds = 5;
Folds = [1:nFolds];
times = 1:length(Y_all);
cvpart = cvpartition(length(Y_all),'kFold',nFolds);                         % create even-ish partitions for k-folding
for iFold = 1:nFolds                                                        % over all Folds
    timesTrain{iFold} = times(cvpart.training(iFold));
    timesTest{iFold}  = times(cvpart.test(iFold));
    Y_train{iFold}    = Y_all(timesTrain{iFold},:);
    Y_test{iFold}     = Y_all(timesTest{iFold},:);
    M_train{iFold}    = Q_all(timesTrain{iFold},:); %M_all(timesTrain{iFold},:);
    M_test{iFold}     = Q_all(timesTest{iFold},:); %M_all(timesTest{iFold},:);
%     M_g{iFold}        = Q_group(timesTrain{iFold},:); %M_group(timesTrain{iFold},:);
%     M_gtest{iFold}    = Q_group(timesTest{iFold},:); %M_group(timesTest{iFold},:);
    nData(iFold)      = cvpart.TestSize(iFold);
end
%% Regularised LS (joint) + CV across all datasets
% addpath('LASSO_Shmidt')                                                     % for lasso library
log_max = 1; log_min = -6; nLambdas = 400;
lambdas = [10.^(log_min + (log_max-log_min)*rand(nLambdas,1)) rand(nLambdas,1)];
fid=fopen('progress.txt','w');
fprintf(fid,'%10s %12s %12s %12s \r\n','RS iter','Gamma','Alpha','RSS');
writeFormat = '%10i %12.4f %12.4f %12.4f \r\n';
for iLambda = 1:nLambdas                                                    % across regiularisation coeffs  
    gamma   = lambdas(iLambda,1);                                           % RLS coefficient (from 10^-6 to 1)
    alpha   = lambdas(iLambda,2);                                           % spasity calibration coefficient (from 0 to 1) 
    for iFold = 1:nFolds                                                    % across folds 
        R_mm  = M_train{iFold}'*M_train{iFold};
        gain     = pinv(R_mm + gamma*eye(size(R_mm)))*M_train{iFold}';      % RLS gain
        g_bar    = gain*Y_train{iFold};                                     % Tikhonov
        g_lasso  = LassoShooting(M_train{iFold},Y_train{iFold},gamma,'verbose',0); % LASSO
        g_spl    = SPLAsso(Y_train{iFold}, M_train{iFold}, p_sparesgroup, (1-alpha)*gamma, alpha*gamma); % sparse group lasso
%% validation for argumented problems
        PE        = Y_test{iFold} - M_test{iFold}*g_bar;
        PE_lasso  = Y_test{iFold} - M_test{iFold}*g_lasso;
        PE_spl    = Y_test{iFold} - M_test{iFold}*g_spl;
        RSS(iLambda,iFold)       = PE'*PE/nData(iFold);
        RSS_lasso(iLambda,iFold) = PE_lasso'*PE_lasso/nData(iFold);
        RSS_spl(iLambda,iFold)  = PE_spl'*PE_spl/nData(iFold);
        clear PE PE_lasso PE_spl       
    end
%% Compute precicted errors and other criteria across all folds
        PRESS(iLambda)          = sum(RSS(iLambda,:));                          
        PRESS_lasso(iLambda)    = sum(RSS_lasso(iLambda,:));
        PRESS_spl(iLambda)      = sum(RSS_spl(iLambda,:));
%% Update progress log
        fid=fopen('progress.txt','a+');
        temp =  [iLambda; gamma; alpha; PRESS_spl(iLambda)];
        fprintf(fid, writeFormat,temp);
        fclose(fid);
end
% find optimal regularisation coefficients
[PRESS_min,i_min] = min(PRESS);
[PRESS_min_l,i_min_l] = min(PRESS_lasso);
[min_val_spl,i_min_spl]=min(PRESS_spl(:));
vec     = [0:1/20:1];
xq      = 10.^(log_min + (log_max-log_min)*vec);
yq      = [0:0.02:1]';
[Xq,Yq] = meshgrid(xq,yq);                                                  % create meshgrid
Z_spgl  = griddata(lambdas(:,1),lambdas(:,2),PRESS_spl,Xq,Yq);              % interpolate function from points
%% Plots
fig('PRESS spgl',visFlag);
mesh(Xq,Yq,Z_spgl);alpha(0.4);
 hold on;
colormap(my_map)
scatter3(lambdas(:,1),lambdas(:,2),PRESS_spl,'filled');
plot3(lambdas(i_min_spl,1),lambdas(i_min_spl,2),min_val_spl,'*','Linewidth',5);
set(gca,'XScale','log')
xlabel('$\gamma$');ylabel('$\alpha$');zlabel('PRESS');
legend('Interpolated PRESS','RS points','Minimum')
title('Sparse group lasso');
tikzName = [folderName,'/PRESS_spgl.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '6cm', 'width','6cm','checkForUpdates',false);
%%        
fig('PRESS',visFlag);
subplot(1,2,1);
scatter(lambdas(:,1),PRESS,'filled'); hold on;
plot(lambdas(i_min,1), PRESS_min, '*','linewidth',3);
set(gca,'XScale','log')
xlabel('$\lambda$');ylabel('PRESS');
title('Tikhonov regularisation');
subplot(1,2,2);
scatter(lambdas(:,1),PRESS_lasso,'filled'); hold on;
plot(lambdas(i_min_l,1), PRESS_min_l, '*','linewidth',5);
set(gca,'XScale','log')
xlabel('$\lambda$');ylabel('PRESS');
title('LASSO regularisation');
xlabel('$\lambda$');ylabel('PRESS');
tikzName = [folderName,'/PRESS_min.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '5cm', 'width','12cm','checkForUpdates',false);
%% Estimating with optimal regularisation coefficient
g_bar   = Q_all\Y_all;
B_bar   = linsolve(R_all,g_bar,struct('UT', true)); 
L       = size(A,2);
Betas_nonreg_opt  = reshape(B_bar,[finalTerm,L]);
%% simple regularisation
lambda_opt = lambdas(i_min,1);
lambda_lasso_opt = lambdas(i_min_l,1);
R_mm    = Q_all'*Q_all;
gain    = inv(R_mm + lambda_opt*eye(size(R_mm)))*Q_all';                    % RLS gain
g_tikh  = gain*Y_all;
B_tikh   = linsolve(R_all,g_tikh,struct('UT', true)); 
g_lasso =  LassoShooting(Q_all,Y_all,lambda_lasso_opt,'verbose',0);
B_lasso   = linsolve(R_all,g_lasso,struct('UT', true)); 
Betas_tikh_opt  = reshape(B_tikh,[finalTerm,L]);
Betas_lasso_opt = reshape(B_lasso,[finalTerm,L]);
%% sparse group lasso
lambda_spl_opt = lambdas(i_min_spl,1);
alpha_spl_opt = lambdas(i_min_spl,2);
g_spl  = SPLAsso(Y_all, Q_all, p_sparesgroup, (1-alpha_spl_opt)*lambda_spl_opt, alpha_spl_opt*lambda_spl_opt); 
B_spl   = linsolve(R_all,g_spl,struct('UT', true)); 
Betas_spl_opt = reshape(B_spl,[finalTerm,L]);
%% display estimates
if disFlag
    Betas_nonreg_opt
    Betas_tikh_opt
    Betas_lasso_opt
    Betas_spl_opt
end
%% Saving external parameters to table
clear Tab
Tab = table(Step,Terms);
for iBeta=1:L
    Parameters = round(Betas_nonreg_opt(:,iBeta),6);
    varName = ['$\beta_{',num2str(iBeta-1),'}$'];
    Tab = addvars(Tab,Parameters,'NewVariableNames',varName);
end
tableName = [folderName,'/Betas_ols'];
table2latex(Tab,tableName);
clear Tab
Tab = table(Step,Terms);
for iBeta=1:L
    Parameters = round(Betas_tikh_opt(:,iBeta),6);
    varName = ['$\beta_{',num2str(iBeta-1),'}$'];
    Tab = addvars(Tab,Parameters,'NewVariableNames',varName);
end
tableName = [folderName,'/Betas_tikhonov'];
table2latex(Tab,tableName);
clear Tab
Tab = table(Step,Terms);
for iBeta=1:L
    Parameters = round(Betas_lasso_opt(:,iBeta),6);
    varName = ['$\beta_{',num2str(iBeta-1),'}$'];
    Tab = addvars(Tab,Parameters,'NewVariableNames',varName);
end
tableName = [folderName,'/Betas_lasso'];
table2latex(Tab,tableName);
clear Tab
Tab = table(Step,Terms);
for iBeta=1:L
    Parameters = round(Betas_spl_opt(:,iBeta),8);
    varName = ['$\beta_{',num2str(iBeta-1),'}$'];
    Tab = addvars(Tab,Parameters,'NewVariableNames',varName);
end
tableName = [folderName,'/Betas_sparse_group'];
table2latex(Tab,tableName);
%% For validation plots
if length(testFiles) > 1
    L1 = 2;
else
    L1 = 1;
end
L2 = round(length(testFiles)/L1);
index_test  = 1:3000;
index_plot  = 1:length(index_test);
%% Validate unconstrained
Theta_test  = Betas_nonreg_opt*A_valid';
iTheta      = 0;
fig('Validation OLS',visFlag);
for iFile = testFiles
% OSA prediction
    fName   = [dictFolder,'/Dict_',dataset,num2str(iFile)];
    File    = matfile(fName,'Writable',true);
    indSign = S(1:finalTerm);                                               % select the indeces of significant terms from the ordered set
    Phi_all = File.term(index_test,:);                                      % extract all terms into a vector
    Phi     = Phi_all(:,indSign);                                           % select only signficant terms
    iTheta  = iTheta + 1;
    y_osa   = Phi*Theta_test(:,iTheta);                                      % model NARMAX output
    RMSE_osa(iTheta) = sqrt(mean((File.y_narx(index_test,1) - y_osa).^2));    % Root Mean Squared Error
% MPO prediction
    y_mpo = File.y_narx(index_test,1);
    u_mpo = File.u_narx(index_test,1);
    for t=n_y+1:index_test(end)
        sum = 0;
        for iTerm=1:finalTerm
            x_mpo(iTerm) = f_mpo{iTerm}(u_mpo,y_mpo,t);
        end
        y_mpo(t)    = x_mpo*Theta_test(:,iTheta);
    end
    RMSE_mpo(iTheta) = sqrt(mean((File.y_narx(index_test,1) - y_mpo).^2)); 
% Compare outputs
    subplot(L2,L1,iTheta);
    plot(index_test(index_plot)+File.t_0,File.y_narx(index_test(index_plot),1),'LineWidth',2); hold on;
    plot(index_test(index_plot)+File.t_0,y_osa(index_plot,1),'--','LineWidth',2); hold on;
    plot(index_test(index_plot)+File.t_0,y_mpo(index_plot,1),'.','LineWidth',3); hold on;
    legend('True output','OSA predition','MPO prediction');
    xlabel('Sample index'); ylabel(['$',y_str,'$']);
    title([dataset,num2str(iFile),': RMSE(OSA) = ',num2str(RMSE_osa(iTheta)),', RMSE(MPO) = ',num2str(RMSE_mpo(iTheta))]);
 clear File Phi_all Phi y_osa y_mpo
end
tikzName = [folderName,'/','OLS_validation.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '12cm', 'width','15cm','checkForUpdates',false);
%% Validate lasso reg
Theta_test  = Betas_lasso_opt*A_valid';
iTheta      = 0;
fig('Validation LAS',visFlag);
for iFile = testFiles
% OSA prediction
    fName   = [dictFolder,'/Dict_',dataset,num2str(iFile)];
    File    = matfile(fName,'Writable',true);
    indSign = S(1:finalTerm);                                               % select the indeces of significant terms from the ordered set
    Phi_all = File.term(index_test,:);                                      % extract all terms into a vector
    Phi     = Phi_all(:,indSign);                                           % select only signficant terms
    iTheta  = iTheta + 1;
    y_osa   = Phi*Theta_test(:,iTheta);                                      % model NARMAX output
    RMSE_osa(iTheta) = sqrt(mean((File.y_narx(index_test,1) - y_osa).^2));    % Root Mean Squared Error
% MPO prediction
    y_mpo = File.y_narx(index_test,1);
    u_mpo = File.u_narx(index_test,1);
    for t=n_y+1:index_test(end)
        sum = 0;
        for iTerm=1:finalTerm
            x_mpo(iTerm) = f_mpo{iTerm}(u_mpo,y_mpo,t);
        end
        y_mpo(t)    = x_mpo*Theta_test(:,iTheta);
    end
    RMSE_mpo(iTheta) = sqrt(mean((File.y_narx(index_test,1) - y_mpo).^2)); 
% Compare outputs
    subplot(L2,L1,iTheta);
    plot(index_test(index_plot)+File.t_0,File.y_narx(index_test(index_plot),1),'LineWidth',2); hold on;
    plot(index_test(index_plot)+File.t_0,y_osa(index_plot,1),'--','LineWidth',2); hold on;
    plot(index_test(index_plot)+File.t_0,y_mpo(index_plot,1),'.','LineWidth',3); hold on;
    legend('True output','OSA predition','MPO prediction');
    xlabel('Sample index'); ylabel(['$',y_str,'$']);
    title([dataset,num2str(iFile),': RMSE(OSA) = ',num2str(RMSE_osa(iTheta)),', RMSE(MPO) = ',num2str(RMSE_mpo(iTheta))]);
      
 clear File Phi_all Phi y_osa y_mpo u_mpo
end
tikzName = [folderName,'/','LASSO_validation.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '12cm', 'width','15cm','checkForUpdates',false);
%% Validate sparse group lasso
% Theta_test  = Betas_spl_opt*A_valid'; % Betas_nonortspl_opt*A_valid'; %
% iTheta      = 0;
% fig('Validation SPGL',visFlag);
% for iFile = testFiles
% % OSA prediction
%     fName   = [dictFolder,'/Dict_',dataset,num2str(iFile)];
%     File    = matfile(fName,'Writable',true);
%     indSign = S(1:finalTerm);                                               % select the indeces of significant terms from the ordered set
%     Phi_all = File.term(index_test,:);                                      % extract all terms into a vector
%     Phi     = Phi_all(:,indSign);                                           % select only signficant terms
%     iTheta  = iTheta + 1;
%     y_osa   = Phi*Theta_test(:,iTheta);                                      % model NARMAX output
%     RMSE_osa(iTheta) = sqrt(mean((File.y_narx(index_test,1) - y_osa).^2));    % Root Mean Squared Error
% % MPO prediction
%     y_mpo = File.y_narx(index_test,1);
%     u_mpo = File.u_narx(index_test,1);
%     for t=n_y+1:index_test(end)
%         sum = 0;
%         for iTerm=1:finalTerm
%             x_mpo(iTerm) = f_mpo{iTerm}(u_mpo,y_mpo,t);
%         end
%         y_mpo(t)    = x_mpo*Theta_test(:,iTheta);
%     end
%     RMSE_mpo(iTheta) = sqrt(mean((File.y_narx(index_test,1) - y_mpo).^2)); 
% % Compare outputs
%     subplot(L2,L1,iTheta);
%     plot(index_test(index_plot)+File.t_0,File.y_narx(index_test(index_plot),1),'LineWidth',2); hold on;
%     plot(index_test(index_plot)+File.t_0,y_osa(index_plot,1),'--','LineWidth',2); hold on;
%     plot(index_test(index_plot)+File.t_0,y_mpo(index_plot,1),'.','LineWidth',3); hold on;
%     legend('True output','OSA predition','MPO prediction');
%     xlabel('Sample index'); ylabel(['$',y_str,'$']);
%     title([dataset,num2str(iFile),': RMSE(OSA) = ',num2str(RMSE_osa(iTheta)),', RMSE(MPO) = ',num2str(RMSE_mpo(iTheta))]);
%  clear File Phi_all Phi y_model
% end
% tikzName = [folderName,'/','Sparse_lasso_validation.tikz'];
% cleanfigure;
% matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
%             false, 'height', '12cm', 'width','15cm','checkForUpdates',false);
%% Validate Tikhonov reg
Theta_test  = Betas_tikh_opt*A_valid';
iTheta = 0;
fig('Validation Tikhonov',visFlag);
for iFile = testFiles
% OSA prediction
    fName   = [dictFolder,'/Dict_',dataset,num2str(iFile)];
    File    = matfile(fName,'Writable',true);
    indSign = S(1:finalTerm);                                               % select the indeces of significant terms from the ordered set
    Phi_all = File.term(index_test,:);                                      % extract all terms into a vector
    Phi     = Phi_all(:,indSign);                                           % select only signficant terms
    iTheta  = iTheta + 1;
    y_osa   = Phi*Theta_test(:,iTheta);                                      % model NARMAX output
    RMSE_osa(iTheta) = sqrt(mean((File.y_narx(index_test,1) - y_osa).^2));    % Root Mean Squared Error
% MPO prediction
    y_mpo = File.y_narx(index_test,1);
    u_mpo = File.u_narx(index_test,1);
    for t=n_y+1:index_test(end)
        sum = 0;
        for iTerm=1:finalTerm
            x_mpo(iTerm) = f_mpo{iTerm}(u_mpo,y_mpo,t);
        end
        y_mpo(t)    = x_mpo*Theta_test(:,iTheta);
    end
    RMSE_mpo(iTheta) = sqrt(mean((File.y_narx(index_test,1) - y_mpo).^2)); 
% Compare outputs
    subplot(L1,L2,iTheta);
    plot(index_test(index_plot)+File.t_0,File.y_narx(index_test(index_plot),1),'LineWidth',2); hold on;
    plot(index_test(index_plot)+File.t_0,y_osa(index_plot,1),'--','LineWidth',2); hold on;
    plot(index_test(index_plot)+File.t_0,y_mpo(index_plot,1),'.','LineWidth',3); hold on;
    legend('True output','OSA predition','MPO prediction');
    xlabel('Sample index'); ylabel(['$',y_str,'$']);
    title([dataset,num2str(iFile),': RMSE(OSA) = ',num2str(RMSE_osa(iTheta)),', RMSE(MPO) = ',num2str(RMSE_mpo(iTheta))]);
 clear File Phi_all Phi y_model
end
tikzName = [folderName,'/','Tikhonov_validation.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '4cm', 'width','15cm','checkForUpdates',false);

%% Save all workspace
foName = ['../SYSDYMATS_data/results/',folderName]
make_folder(foName)
fiName = [foName,'/','Contrained_estimation_resullts'];
save(fiName);

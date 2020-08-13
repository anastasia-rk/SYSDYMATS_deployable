clear all; local_init;
%% Select the dataset and domain
dataset = 'V';  
folder = 'results_cv';                                                     % specify category where to save files
dFolder = 'dictionaries';
regressors = questdlg('Select the domain of regressors', ...
    'Domain choice',...
	'Shift','Delta','');
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
        metaFileName = ['Meta_delta_',dataset];
        load(metaFileName);
         folder = ['delta_',folder];
         dFolder = ['delta_',dFolder];
         if normC ~= 1
            folder = [folder,'_norm'];
            dFolder = [dFolder,'_norm'];
         end
        direction = questdlg('Type of delta operator', ...
        'Causality',...
        'Forward','Backward','');
        switch direction
            case 'Backward'
                folder = [folder,'_b'];
                dFolder = [dFolder,'_b'];
            case 'Forward'
                folder = [folder,'_f'];
                dFolder = [dFolder,'_f'];
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
Files_all =  1:K;                                                           % ids of the sample files
Files = Files_all(2:end-1);
testFiles   = [1 4 6 K];
K = length(Files);
% Set maximum number of covariates
if length(dict_terms) < 30                                                  % Maximum significant terms (if the algorithm is not terminated by the criterion)
    maxSign = length(dict_terms);
else
    maxSign = 30;                                                           
end
dict_terms_all = dict_terms;
%% Create the parameter matrix based on the dataset
load(['External_parameters_',dataset]);
x = values(Files,1);                                                        % design parameter values for training
x_valid = values(testFiles,1);                                              % design parameter values for validation
% x_m = mean(x);
% x_s = sqrt((x-x_m)'*(x-x_m));
% x   = x/x_s;
if size(values,2) > 1
    y = values(Files,2);
    y_valid = values(testFiles,2);
%     y_m = mean(y);
%     y_s = sqrt((y-y_m)'*(y-y_m));
%     y   = y/y_s;
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
%% Preliminaries for structure identification
index       = 1:nNarx;                                                      % structure id over all times
dict_terms  = dict_terms_all;                                               % refill the dictionary
%% Select first significant basis vector for all datasets
iTerm = 1;                                                                  % the first significant term
AERR{iTerm} = zeros(nTerms,1);                                              % placeholder for AERR criteria
for iFile = Files                                                           % over all datasets
    fName   = [dictFolder,'/',char(fileNames(iFile))];
    File    = matfile(fName,'Writable',true);
    y_n     = File.y_narx(:,1);
    residual_init{iFile} =  y_n(index,1);                                   % initial residual
    for jTerm = dict_terms                                                  % over all polynomial terms in the dictionary
        term_all    = File.term(:,jTerm);
        term0       = term_all(index,:);
        cf(iFile,jTerm)     = cor_sqr(residual_init{iFile},term0);          % squared correlation coefficient for the dataset and the polynomial term
        AERR{iTerm}(jTerm)  = AERR{iTerm}(jTerm) + cf(iFile,jTerm);         % Average error reduction ration over all datasets
        clear term0 term_all
    end
    clear File y_n 
end
AERR{iTerm}(:,:)    = AERR{iTerm}(:,:)/K;
[AERR_m,iMax]       = max(AERR{iTerm});                                     % find the index of the term with the highest criterion across all datasets
AERR_mm(iTerm,1)    = AERR_m;
S(iTerm)            = iMax;                                                 % save index of the term
dict_terms(iMax)    = [];                                                   % reduce the dictionary of available term
BIC_sum             = 0;
for iFile = Files                                                           % over all datasets
    fName = [dictFolder,'/',char(fileNames(iFile))];
    File  = matfile(fName,'Writable',true);
    term_all    = File.term(:,iMax);
    alpha{iFile}(:,iTerm)    = term_all(index,:);                           % the corresponding basis candidate term    
    phi  {iFile}(:,iTerm)    = term_all(index,:);                           % the corresponding basis vector 
    residual{iFile}(:,iTerm) = residual_update(residual_init{iFile},...     % the corresponding model residual
                                               phi{iFile}(:,iTerm));                                                        
    BIC_sum  = BIC_sum  +  BIC(residual{iFile}(:,iTerm),nNarx,iTerm);       % BIC for the iFile dataset
    clear File term_all
end
    BIC_all(iTerm)            = BIC_sum/K;                                  % average AMDL over all sets
    significant_term{iTerm}   = symb_term{S(iTerm)};
%% Identify optimal model structure - EFOR loop  
converged   = false;
bics        = [];
while(iTerm < maxSign) %&& ~converged                                       % loop over the number of significant terms
     iTerm = iTerm + 1;                                                     % increase the number of significant terms
     AERR{iTerm} = zeros(nTerms,1);                                         % placeholder for AERR criteria
     for iFile = Files                                                      % over all datasets
         fName   = [dictFolder,'/',char(fileNames(iFile))];
         File    = matfile(fName,'Writable',true);
         for jTerm = dict_terms                                             % over all polynomial terms in the dictionary
             term_all    = File.term(:,jTerm);
             p{iTerm,iFile}(:,jTerm) = orthogonalise(term_all(index,:),...
                                                    phi{iFile},iTerm);      % orthogonalise basis
             cf(iFile,jTerm)         = cor_sqr(residual_init{iFile},...
                                              p{iTerm,iFile}(:,jTerm));     % squared correlation coefficient for the dataset and the polynomial term
             AERR{iTerm}(jTerm) = AERR{iTerm}(jTerm) + cf(iFile,jTerm);     % average error reduction ration over all datasets
         end
         clear File
     end
     AERR{iTerm}(:,:)    = AERR{iTerm}(:,:)/K;
     [AERR_m,iMax]       = max(AERR{iTerm});                                % find the index of the term with the highest criterion across all datasets
     AERR_mm(iTerm,1)    = AERR_m;
     S(iTerm)            = iMax;                                            % save index of the term  
     ind = find(dict_terms == iMax);
     dict_terms(ind) = [];                                                  % Reduce the dictionary of available terms
     BIC_sum         = 0;
     for iFile = Files
         fName   = [dictFolder,'/',char(fileNames(iFile))];
         File    = matfile(fName,'Writable',true);
         alpha{iFile}(:,iTerm) = File.term(index,S(iTerm));                 % the corresponding basis candidate term    
         phi{iFile}(:,iTerm)   = p{iTerm,iFile}(:,S(iTerm));                % the corresponding basis vector 
         residual{iFile}(:,iTerm) = residual_update(residual{iFile}(:,iTerm-1),...
                                                   phi{iFile}(:,iTerm));    % the corresponding model residual                                 
         BIC_sum  = BIC_sum  +  BIC(residual{iFile}(:,iTerm),nNarx,iTerm);  % BIC for the iFile dataset
         clear File x_n
     end
     significant_term{iTerm} = symb_term{S(iTerm)};
     BIC_all(iTerm) = BIC_sum/K;                                            % average AMDL over all sets
     converged_BIC  = (abs((BIC_all(iTerm) - BIC_all(iTerm-1))/BIC_all(iTerm)) < 0.002); % check convergence
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
%% Estimate internal polynomial coefficients
for iFile=Files
    U{iFile} = zeros(finalTerm,finalTerm);                                  % placeholder for upper-trig unit matrix
    iTerm = 1;                                                              % for the first term
    g{iFile}(iTerm) = (residual_init{iFile}'*phi{iFile}(:,iTerm))/...
                      (phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));           % top raw of the rotation matrix
    for jTerm =iTerm:finalTerm
        U{iFile}(iTerm,jTerm) = alpha{iFile}(:,jTerm)'*phi{iFile}(:,iTerm)/...
                                (phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));
    end
    for iTerm = 2:finalTerm                                                 % loop over significant terms
        g{iFile}(iTerm,1) = (residual{iFile}(:,iTerm-1)'*phi{iFile}...
                    (:,iTerm))/(phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));  % righthand side (normalised)
        for jTerm =iTerm:finalTerm
            U{iFile}(iTerm,jTerm) = alpha{iFile}(:,jTerm)'*phi{iFile}...
                     (:,iTerm)/(phi{iFile}(:,iTerm)'*phi{iFile}(:,iTerm));  % upper triangular unit matrix
        end
    end
    Dets = det(U{iFile});
    Theta(:,iFile) = linsolve(U{iFile},g{iFile},struct('UT', true));        % solve upper triangular system via backward substitution
    Parameters = round(Theta(:,iFile),2);
    varName = [dataset,num2str(iFile)];
    Tab = addvars(Tab,Parameters,'NewVariableNames',varName);
end
%% Save internal parameters to table
Tab             = addvars(Tab,AERR,'NewVariableNames',{'AERR($\%$)'});
internalParams  = addvars(Tab,BIC_trunc','NewVariableNames',{'BIC'});
if disFlag 
   internalParams
end
tableName = [folderName,'/Thetas_overall'];
table2latex(internalParams,tableName);
clear Tab tableName
clear AERR alpha U phi residual p g
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
splits  = 2;                                                                 % split time between training and testing data                       
iMpo   = 0;
for iFile = Files
    fName     = [dictFolder,'/',char(fileNames(iFile))];
    File      = matfile(fName,'Writable',true);
    indSign   = S(1:finalTerm);                                             % select the indeces of significant terms from the ordered set
    Phi_file  = File.term(:,:);                                             % extract all terms into a vector - cannot reorder directly in files
    y_file    = File.y_narx(:,:);                                           % extract output
    Phi       = Phi_file(:,indSign);                                        % select only signficant terms
    Phi_bar   = blkdiag(Phi_bar,Phi);                                       % diagonal block, size TKxNK
    Y_all     = [Y_all; y_file];                                            % block vector, size TKx1
    iMpo      = iMpo + 1;
    t_full(iMpo+1) = length(y_file);
    t_test{iMpo}   = 1:length(y_file); %[(splits-1)*floor(length(y_file)/splits):length(y_file)];
end
M_all  = Phi_bar*Kr;  
[Q_all,R_all] = mgson(M_all);                                               % Orthogonalise M via modified Gram-Schmidt
% !!!!!!!!!!!! From this point onwards we work with Q_all for
% cross-validation. We will utilise the R_all for orthogonal solution only after
% optimal regularisation parameters have been found
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
splits  = 2;                                                                 % split time between training and testing data                       
iMpo   = 0;
for iFile = Files
    fName     = [dictFolder,'/',char(fileNames(iFile))];
    File      = matfile(fName,'Writable',true);
    indSign   = S(1:finalTerm);                                             % select the indeces of significant terms from the ordered set
    Phi_file  = File.term(:,:);                                             % extract all terms into a vector - cannot reorder directly in files
    y_file    = File.y_narx(:,:);                                           % extract output
    Phi       = Phi_file(:,indSign);                                        % select only signficant terms
    Phi_bar   = blkdiag(Phi_bar,Phi);                                       % diagonal block, size TKxNK
    Y_all     = [Y_all; y_file];                                            % block vector, size TKx1
    iMpo      = iMpo + 1;
    t_full(iMpo+1) = length(y_file);
    t_test{iMpo}   = 1:length(y_file); %[(splits-1)*floor(length(y_file)/splits):length(y_file)];
end
M_all  = Phi_bar*Kr;  
[Q_all,R_all] = mgson(M_all);                                               % Orthogonalise M via modified Gram-Schmidt
% !!!!!!!!!!!! From this point onwards we work with Q_all for
% cross-validation. We will utilise the R_all for orthogonal solution only after
% optimal regularisation parameters have been found
%% Form folds for CV with K-folds
nFolds = length(Files);
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
    nData(iFold)      = cvpart.TestSize(iFold);
end
%% Form blocks for CV of time series
fig('MPO CV folds',visFlag);
for iMpo = 1:length(Files)
    t_train = times;
    index   = sum(t_full(1:iMpo)) + t_test{iMpo};                           % index of test dataset in the overall time sequence
    t_train(index) = [];                                                    % remove the points of the testing dataset
    iFile   = Files(iMpo);
    Q_train_mpo{iFile}  = M_all(t_train,:);
    Y_train_mpo{iFile}  = Y_all(t_train,:); 
    nData_mpo(iFile)    = length(t_test{iMpo});
    plot(t_train,iMpo*ones(size(t_train)),'.g','Linewidth',5); hold on;
    plot(index,iMpo*ones(size(t_test{iMpo})),'.k','Linewidth',5); hold on;
    clear t_train
end
xlabel('Indeces'); ylabel('Blocks')
%% Random search vector
addpath('LASSO_Shmidt')                                                     % for lasso library
log_max = 0; log_min = -5; nLambdas = 50;
lambdas = [10.^(log_min + (log_max-log_min)*rand(nLambdas,1)) rand(nLambdas,1)];
L       = size(A,2);
fid_osa=fopen([dataset,'_progress_OSA.txt'],'w');
fprintf(fid_osa,'%10s %12s %12s %12s %12s \r\n','RS iter','Gamma','Alpha','PRESS','T elapsed');
fid_mpo=fopen([dataset,'_progress_MPO.txt'],'w');
fprintf(fid_mpo,'%10s %12s %12s %12s %12s \r\n','RS iter','Gamma','Alpha','PRESS','T elapsed');
writeFormat = '%10i %12.4f %12.4f %12.4f %12.4f\r\n';
%% OSA cross-validation for constrained problems
if disFlag
   fprintf(['OSA CV... \n'])
end
pool =  parpool('local');
for iLambda = 1:nLambdas                                                    % across regiularisation coeffs  
    lambda   = lambdas(iLambda,1);                                          % RLS coefficient (from 10^-6 to 1)
    lambda_g = lambdas(iLambda,2);                                          % spasity calibration coefficient (from 0 to 1) 
    tic;
    parfor iFold = 1:nFolds                                                 % across folds 
        R_mm  = M_train{iFold}'*M_train{iFold};
        gain     = inv(R_mm + lambda*eye(size(R_mm)))*M_train{iFold}';      % RLS gain
        g_bar    = gain*Y_train{iFold};                                     % Tikhonov
        g_lasso  = LassoShooting(M_train{iFold},Y_train{iFold},lambda,'verbose',0); % LASSO
%         g_spl    = SPLAsso(Y_train{iFold}, M_train{iFold}, p_sparesgroup, (1-lambda_g)*lambda, lambda_g*lambda); % sparse group lasso
% Validation
        PE        = Y_test{iFold} - M_test{iFold}*g_bar;
        PE_lasso  = Y_test{iFold} - M_test{iFold}*g_lasso;
%         PE_spl    = Y_test{iFold} - M_test{iFold}*g_spl;
        RSS(iLambda,iFold)        = PE'*PE/nData(iFold);
        RSS_lasso(iLambda,iFold)  = PE_lasso'*PE_lasso/nData(iFold);
%         RSS_spl(iLambda,iFold)    = PE_spl'*PE_spl/nData(iFold);
    end
    t_el_osa = toc;
    PRESS(iLambda)       = sum(RSS(iLambda,:));                          
    PRESS_lasso(iLambda) = sum(RSS_lasso(iLambda,:));
%     PRESS_spl(iLambda)   = sum(RSS_spl(iLambda,:));
%% Update progress log
        fid=fopen([dataset,'_progress_OSA.txt'],'a+');
        temp =  [iLambda; lambda; lambda_g; PRESS_lasso(iLambda); t_el_osa];
        fprintf(fid, writeFormat,temp);
        fclose(fid);
        clear temp fid
    if disFlag
        fprintf([num2str(iLambda) ' done... \n'])
    end
end
%% MPO cross-validation for constrained problems
if disFlag
   fprintf(['MPO CV... \n'])
end
for iLambda = 1:nLambdas                                                    % across regiularisation coeffs  
    lambda   = lambdas(iLambda,1);                                          % RLS coefficient (from 10^-6 to 1)
    lambda_g = lambdas(iLambda,2);                                          % spasity calibration coefficient (from 0 to 1) 
    tic;
    parfor  iMpo =1:length(Files) % for MPO - instead of folds loop over files. 
        iFile    = Files(iMpo); % counter to go through design parameter matrix
        R_mm     = Q_train_mpo{iFile}'*Q_train_mpo{iFile};
        gain     = inv(R_mm + lambda*eye(size(R_mm)))*Q_train_mpo{iFile}';         % RLS gain
        g_bar    = gain*Y_train_mpo{iFile};                                        % Tikhonov
        g_lasso  = LassoShooting(Q_train_mpo{iFile},Y_train_mpo{iFile},lambda,'verbose',0); % LASSO
        g_spl    = SPLAsso(Y_train_mpo{iFile}, Q_train_mpo{iFile}, p_sparesgroup, (1-lambda_g)*lambda, lambda_g*lambda); % sparse group lasso
        Betas_tikh   = reshape(g_bar); % reshape(linsolve(R_all,g_bar,struct('UT', true)),[finalTerm L]); 
        Betas_lass   = reshape(g_lasso); %reshape(linsolve(R_all,g_lasso,struct('UT', true)),[finalTerm L]);  
        Betas_spgl   = reshape(g_spl); %reshape(linsolve(R_all,g_spl,struct('UT', true)),[finalTerm L]);  
        Theta_tikh{iMpo}  = Betas_tikh*A(iMpo,:)';
        Theta_lass{iMpo}  = Betas_lass*A(iMpo,:)';
        Theta_spgl{iMpo}  = Betas_spgl*A(iMpo,:)';
    end
    t_el_mpo = toc;
% MPO validation
    iMpo = 0;
    for iFiles = Files
        iMpo    = iMpo + 1;
        fName   = [dictFolder,'/dict_',dataset,num2str(iFile)];
        File    = matfile(fName,'Writable',true);
        y_mpo   = File.y_narx(t_test{iMpo},1);
        y_mpo_lasso = y_mpo;
        y_mpo_spgl  = y_mpo;
        u_mpo = File.u_narx(t_test{iMpo},1);
        for t=n_y+1:nData_mpo(iFile)
            for iTerm=1:finalTerm
                x_mpo(iTerm)        = f_mpo{iTerm}(u_mpo,y_mpo,t);
                x_mpo_lasso(iTerm)  = f_mpo{iTerm}(u_mpo,y_mpo_lasso,t);
                x_mpo_spgl(iTerm)   = f_mpo{iTerm}(u_mpo,y_mpo_spgl,t);
            end
            y_mpo(t)        = x_mpo*Theta_tikh{iMpo};
            y_mpo_lasso(t)  = x_mpo_lasso*Theta_lass{iMpo};
            y_mpo_spgl(t)   = x_mpo_spgl*Theta_spgl{iMpo};
        end
        PE_mpo        = File.y_narx(t_test{iMpo},1) - y_mpo;
        PE_lasso_mpo  = File.y_narx(t_test{iMpo},1) - y_mpo_lasso;
        PE_spl_mpo    = File.y_narx(t_test{iMpo},1) - y_mpo_spgl;
        RSS_mpo(iLambda,iMpo)        = PE_mpo'*PE_mpo/nData_mpo(iFile);
        RSS_lasso_mpo(iLambda,iMpo)  = PE_lasso_mpo'*PE_lasso_mpo/nData_mpo(iFile);
        RSS_spl_mpo(iLambda,iMpo)    = PE_spl_mpo'*PE_spl_mpo/nData_mpo(iFile);
    end
    PRESS_mpo(iLambda)       = sum(RSS_mpo(iLambda,:));                          
    PRESS_lasso_mpo(iLambda) = sum(RSS_lasso_mpo(iLambda,:));
    PRESS_spl_mpo(iLambda)   = sum(RSS_spl_mpo(iLambda,:));
%% Update progress log
        fid=fopen([dataset,'_progress_MPO.txt'],'a+');
        temp =  [iLambda; lambda; lambda_g; PRESS_spl_mpo(iLambda); t_el_mpo];
        fprintf(fid, writeFormat,temp);
        fclose(fid);
        clear temp fid
    if disFlag
        fprintf([num2str(iLambda) ' done... \n'])
    end
end
delete(pool);
clear M_train Y_train Q_train_mpo Y_train_mpo times M_test Y_test 
%% Find optimal regularisation coefficients - CV loop
[PRESS_min,i_min]       = min(PRESS);
[PRESS_min_l,i_min_l]   = min(PRESS_lasso);
% [min_val_spl,i_min_spl] = min(PRESS_spl);
[PRESS_min_mpo,i_min_mpo]       = min(PRESS_mpo);
[PRESS_min_l_mpo,i_min_l_mpo]   = min(PRESS_lasso_mpo);
[min_val_spl_mpo,i_min_spl_mpo] = min(PRESS_spl_mpo);
vec     = [0:1/50:1];
xq      = 10.^(log_min + (log_max-log_min)*vec);
yq      = [0:0.02:1]';
[Xq,Yq] = meshgrid(xq,yq);                                  % create meshgrid
% Z_spgl_osa  = griddata(lambdas(:,1),lambdas(:,2),PRESS_spl,Xq,Yq);
Z_spgl_mpo  = griddata(lambdas(:,1),lambdas(:,2),PRESS_spl_mpo,Xq,Yq);
%% Plot spgl PRESS
% fig('PRESS OSA',visFlag);
% mesh(Xq,Yq,Z_spgl_osa);alpha(0.4);
% hold on;
% colormap(my_map)
% scatter3(lambdas(:,1),lambdas(:,2),PRESS_spl,'filled');
% plot3(lambdas(i_min_spl,1),lambdas(i_min_spl,2),min_val_spl,'*','Linewidth',5);
% set(gca,'XScale','log')
% xlabel('$\gamma$');ylabel('$\alpha$');zlabel('PRESS');
% legend('Interpolated PRESS','RS points','Minimum','Location','northwest')
% print('Spgl_press_osa','-depsc')
% tikzName = [folderName,'/PRESS_spgl_osa.tikz'];
% cleanfigure;
% matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
%             false, 'height', '10cm', 'width','10cm','checkForUpdates',false);
fig('PRESS MPO',visFlag);
mesh(Xq,Yq,Z_spgl_mpo);alpha(0.4);
hold on;
colormap(my_map)
scatter3(lambdas(:,1),lambdas(:,2),PRESS_spl_mpo,'filled');
plot3(lambdas(i_min_spl_mpo,1),lambdas(i_min_spl_mpo,2),min_val_spl_mpo,'*','Linewidth',5);
set(gca,'XScale','log')
xlabel('$\gamma$');ylabel('$\alpha$');zlabel('PRESS');
legend('Interpolated PRESS','RS points','Minimum','Location','northwest')
print('Spgl_press_mpo','-depsc')
tikzName = [folderName,'/PRESS_spgl_mpo.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '10cm', 'width','10cm','checkForUpdates',false);

%% Plot ridge and lasso PRESS     
fig('PRESS',visFlag);
subplot(2,2,1);
scatter(lambdas(:,1),PRESS,'filled'); hold on;
plot(lambdas(i_min,1), PRESS_min, '*','linewidth',3);
set(gca,'XScale','log')
legend('RS points','min')
xlabel('$\gamma$');ylabel('PRESS');
title('Tikhonov regularisation OSA');
subplot(2,2,3)
scatter(lambdas(:,1),PRESS_mpo,'filled'); hold on;
plot(lambdas(i_min_mpo,1), PRESS_min_mpo, '*','linewidth',3);
set(gca,'XScale','log')
legend('RS points','min')
xlabel('$\gamma$');ylabel('PRESS');
title('Tikhonov regularisation MPO');
subplot(2,2,2);
scatter(lambdas(:,1),PRESS_lasso,'filled'); hold on;
plot(lambdas(i_min_l,1), PRESS_min_l, '*','linewidth',5);
set(gca,'XScale','log')
legend('RS points','min')
xlabel('$\lambda$');ylabel('PRESS');
title('LASSO regularisation OSA');
xlabel('$\gamma$');ylabel('PRESS');
subplot(2,2,4);
scatter(lambdas(:,1),PRESS_lasso_mpo,'filled'); hold on;
plot(lambdas(i_min_l_mpo,1), PRESS_min_l_mpo, '*','linewidth',3);
set(gca,'XScale','log')
legend('RS points','min')
xlabel('$\lambda$');ylabel('PRESS');
title('LASSO regularisation MPO');
xlabel('$\gamma$');ylabel('PRESS');
tikzName = [folderName,'/PRESS_min.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '5cm', 'width','12cm','checkForUpdates',false);
%% Estimate parameters with optimal regularisation coefficient
g_bar   = Q_all\Y_all;
B_bar   = linsolve(R_all,g_bar,struct('UT', true)); 
Betas_nonreg_opt  = reshape(B_bar,[finalTerm,L]);
%% Simple regularisation
lambda_opt          = lambdas(i_min,1);
lambda_lasso_opt    = lambdas(i_min_l,1);
lambda_opt_mpo      = lambdas(i_min_mpo,1);
lambda_lasso_opt_mpo = lambdas(i_min_l_mpo,1);
R_mm    = Q_all'*Q_all;
gain    = inv(R_mm + lambda_opt_mpo*eye(size(R_mm)))*Q_all';                    % RLS gain
g_tikh  = gain*Y_all;
B_tikh  = linsolve(R_all,g_tikh,struct('UT', true)); 
g_lasso =  LassoShooting(Q_all,Y_all,lambda_lasso_opt_mpo,'verbose',0);
B_lasso = linsolve(R_all,g_lasso,struct('UT', true)); 
Betas_tikh_opt  = reshape(B_tikh,[finalTerm,L]);
Betas_lasso_opt = reshape(B_lasso,[finalTerm,L]);
%% Sparse group lasso
% lambda_spl_opt = lambdas(i_min_spl,1);
% alpha_spl_opt = lambdas(i_min_spl,2);
lambda_spl_mpo = lambdas(i_min_spl_mpo,1);
alpha_spl_mpo = lambdas(i_min_spl_mpo,2);
g_spl  = SPLAsso(Y_all, Q_all, p_sparesgroup, (1-alpha_spl_mpo)*lambda_spl_mpo, alpha_spl_mpo*lambda_spl_mpo); 
B_spl   = linsolve(R_all,g_spl,struct('UT', true)); 
Betas_spl_opt = reshape(B_spl,[finalTerm,L]);
%% display estimates
if disFlag
    Betas_nonreg_opt
    Betas_tikh_opt
    Betas_lasso_opt
    Betas_spl_opt
end
%% Save all workspace
foName = ['../SYSDYMATS_data/results/',folderName]
make_folder(foName)
fiName = [foName,'/','Contrained_estimation_resullts'];
save(fiName);
fiYunpeng  = ['Results_for_Yunpeng_',dataset];
save(fiYunpeng,'Betas_nonreg_opt','Betas_tikh_opt','Betas_lasso_opt','Betas_spl_opt','A','A_valid','Files','testFiles','A_symb','f_mpo','Terms');
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
    fName   = [dictFolder,'/dict_',dataset,num2str(iFile)];
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
    fName   = [dictFolder,'/dict_',dataset,num2str(iFile)];
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
Theta_test  = Betas_spl_opt*A_valid'; % Betas_nonortspl_opt*A_valid'; %
iTheta      = 0;
fig('Validation SPGL',visFlag);
for iFile = testFiles
% OSA prediction
    fName   = [dictFolder,'/dict_',dataset,num2str(iFile)];
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
 clear File Phi_all Phi y_model
end
tikzName = [folderName,'/','Sparse_lasso_validation.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '12cm', 'width','15cm','checkForUpdates',false);
%% Validate Tikhonov reg
Theta_test  = Betas_tikh_opt*A_valid';
iTheta = 0;
fig('Validation Tikhonov',visFlag);
for iFile = testFiles
% OSA prediction
    fName   = [dictFolder,'/dict_',dataset,num2str(iFile)];
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
            false, 'height', '12cm', 'width','15cm','checkForUpdates',false);
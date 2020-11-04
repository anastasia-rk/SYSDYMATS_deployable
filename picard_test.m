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
        dFolder = ['../SYSDYMATS_dictionaries/',dFolder];
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
         dFolder = ['../SYSDYMATS_dictionaries/',dFolder];
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
       A_symb{iCol+1} = ['$x^{',num2str(powers(iCol,1)),'}$ $y^{',num2str(powers(iCol,2)),'}$'];
       temp_g = ['xi(1)^(',num2str(powers(iCol,1)),')*xi(2)^(',num2str(powers(iCol,2)),')'];
       g_model{iCol+1} = inline(temp_g,'xi'); %@(x,power) x(1)^(power(1))*x(2)^(power(2));
       clear temp_g
   end
else                                                                        % unknown mapping is a curve
    powers = [1 2 3 -1 -2 -3];
    nCols =  min(length(powers),K-1);                                       % number of terms in the model shouldn't be higher then K
    for iCol = 1:nCols                                                      % limit order of the model by the number of experimants
       A = [A x.^powers(iCol)];
       A_valid = [A_valid x_valid.^powers(iCol)];
       A_symb{iCol+1} = ['$x^{',num2str(powers(iCol)),'}$'];
       temp_g = ['xi^(',num2str(powers(iCol)),')'];
       g_model{iCol+1} = inline(temp_g,'xi'); % @(x,power) x^(power);
       clear temp_g
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
    f_model{iTerm} = inline(inline_term{S(iTerm)},'u','y','t');
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
%% Save internal parameters to table
Tab             = addvars(Tab,AERR,'NewVariableNames',{'AERR($\%$)'});
internalParams  = addvars(Tab,BIC_trunc','NewVariableNames',{'BIC'});
if disFlag 
   internalParams
end
tableName = [folderName,'/Terms_efor'];
table2latex(internalParams,tableName);
clear Tab tableName
clear AERR alpha U phi residual p g
%% Form regression matrices for all time points
t_span = 1:200;
for iTerm = 1:finalTerm
% vectorisation for joint estimation
    I   = eye(iTerm);                                                       % unit matrix, size NxN
    Kr  = kron(A,I);
    L   = size(A,2);
    Phi_bar = [];
    Y_all   = [];
    splits  = 2;                                                                 % split time between training and testing data                       
    iMpo   = 0;
    for iFile = Files
        fName     = [dictFolder,'/',char(fileNames(iFile))];
        File      = matfile(fName,'Writable',true);
        indSign   = S(1:iTerm);                                             % select the indeces of significant terms from the ordered set
        Phi_file  = File.term(:,:);                                             % extract all terms into a vector - cannot reorder directly in files
        y_file    = File.y_narx(t_span,:);                                           % extract output
        Phi       = Phi_file(t_span,indSign);                                        % select only signficant terms
        Phi_bar   = blkdiag(Phi_bar,Phi);                                       % diagonal block, size TKxNK
        Y_all     = [Y_all; y_file];                                            % block vector, size TKx1
    end
    M = Phi_bar*Kr;
     %% Picard computation
%         [LSV,SIG,RSV] = svd(Phi_bar,'econ');
        [UV,VV,XV,CV,SV] = gsvd(M,eye(size(M,2)),0); % this is for a regularised problem
        gamma =  sqrt(diag(CV'*CV)./diag(SV'*SV));
        numer = abs(UV'*Y_all);
        Gen_sv{iTerm} = flip(gamma);
        Fourtier{iTerm} = flip(numer);
%         denom = diag(SIG);
        Picard{iTerm} = numer./flip(gamma);
        fig([significant_term{iTerm}],'on');
        subplot(2,2,1)
        semilogy(Gen_sv{iTerm},'o'); hold on;
        semilogy(Fourtier{iTerm},'x'); hold on;
        semilogy(Picard{iTerm}); hold on;
        legend('$\gamma_i$','$|u_i^TY|$','$|u_i^TY|/\gamma_i$');
        xlabel('i');ylabel('GSVD decompostion');
        % Picard per Hansen paper
        [UV1,SV1,VV1] = svd(M,'econ');
        numer1 = abs(UV1'*Y_all);
        sigma  = SV1;
        subplot(2,2,2)
        Picard_Hansen{iTerm} = picard(UV1,Gen_sv{iTerm},Y_all,1);
%% Condition numbers - effective etc
        condition_normal{iTerm} = cond(A);
        condition_effective{iTerm} = (1/SV1(end))*norm(numer1)/norm(pinv(diag(Gen_sv{iTerm}))*numer1);
        lse = pinv(M)*Y_all;
        for iComponent=1:size(M,2)
            condition_component{iTerm}(iComponent) = norm(M(:,iComponent))*norm(Y_all)/abs(lse(iComponent));
        end
        subplot(2,2,3)
        semilogy(condition_component{iTerm}); hold on;
        semilogy(1:iComponent,condition_effective{iTerm}*ones(1,iComponent)); hold on;
        semilogy(1:iComponent,condition_normal{iTerm}*ones(1,iComponent)); hold on;
        legend('$\rho(X,Y,b_i)$','$\eta(X,Y)$','$\kappa(X)$');
        xlabel('i');ylabel('Condition numbers');
        index = find()
end
% %% Plots
% fig('Picard gsvd','on')
% for iTerm=1:finalTerm
%     semilogy(Picard{iTerm}); hold on;
% %     legends{iTerm} = char(significant_term{iTerm});
% end
% 
% fig('Picard svd','on')
% for iTerm=1:finalTerm
%     semilogy(Picard_Hansen{iTerm}); hold on;
%     legends{iTerm} = significant_term{iTerm};
% end
%%
[coeffs,FitInfo] = lasso(M,Y_all,'lambda',0.0001);
lassoPlot(coeffs,FitInfo,'PlotType','CV');
legend('show')
figure;
Y_new = M*coeffs;
plot(Y_all); hold on;
plot(Y_new)
hold off

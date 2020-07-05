function[done] = delta_narx(foamset,dataset,direction,normalise)
% Creates dictionaries of delta-domain regressors from the time-series data
local_init;
%% Set up global parameters
difFormat = false;
switch foamset
    case 'foam_2010'
        input_i  = 3;                                                       % input column index
        output_i = 2;                                                       % output column index
        K        = 10;                                                      % number of datasets
        normC    = 400;
        dT       = 0.02;
    case 'foam_2019'
        input_i  = 2;                                                       % input column index
        output_i = 3;                                                       % output column index
        K        = 12;                                                      % number of datasets
        normC    = 100;
        dT       = 2.3438e-04;
        difFormat = true;
    case 'foam_2020'
        input_i  = 2;                                                       % input column index
        output_i = 3;                                                       % output column index
        dT  = 9/51200;                                                      % sampling time
        switch dataset
            case 'HS'
                K = 30; 
            case 'VS'
                K = 30; 
            case 'HF'
                K = 60;                                         
        end
        normC    = 1;
         difFormat = true;
    case 'sweep_2020'
        input_i  = 2;                                                       % input column index
        output_i = 3;                                                       % output column index
        dT       = 9/51200;                                                 % sampling time
        K        = 6;
        normC    = 1;
end
switch normalise
    case 'No'
        folder = 'delta_dictionaries';                                      % specify category where to save files
        normC = 1;
    case 'Yes'
        folder = 'delta_dictionaries_norm';                                 % specify category where to save files
end
switch direction
    case 'Backward'
        folder = [folder,'_b'];
    case 'Forward'
        folder = [folder,'_f'];
end
addpath(foamset)
% Length of input and output lags
lambda  = 2;                                                                % order of delta operator
lambda_p = 3;                                                               % order of polynomial
names = {'set','lambda'};                                                   % names used to define results folder name (no more than 3).
folderName = make_folder(folder,names,dataset,lambda);                      % create results folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create string array of input vector
switch direction
    case 'Backward'
        folder = [folder,'_b'];
        for it=1:lambda
            x_str{it} = ['\delta^',num2str(it),' y(t)'];
        end
         x_str{lambda+1} = 'u(t-1)';
         y_str = 'y(t)';
         poly_index = lambda+1;
         for it=lambda+2:2*lambda
            x_str{it} = ['\delta^',num2str(it-lambda-1),' u(t-1)'];
         end
    case 'Forward'
         x_str{1} = 'y(t)';
         for it=2:lambda
            x_str{it} = ['\delta^',num2str(it-1),' y(t)'];
         end
         x_str{lambda+1} = 'u(t-1)';
         for it=lambda+2:2*lambda
            x_str{it} = ['\delta^',num2str(it-lambda-1),' u(t-1)'];
         end
         poly_index = [1 lambda+1];
         y_str = ['\delta^',num2str(lambda),' y(t)'];
end
d = it;
t_0 = lambda+5;
T   = length(Input);                                                        % length of the observation sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create sum index permutations
% Intial combination accounts for all linear regressors
indeces{1} = [1:d]';
% Products with only unique elements
for iLambda=2:lambda_p                                                      % iLambda - order of polynomial term
    M = permn(1:d,iLambda);                                                 % get all permutations with repetition
    for j=iLambda:-1:2
        ind = find(M(:,j)>M(:,j-1));                                        % sort out only increasing indeces
        M = M(ind,:);                                                       % update set
        clear ind
    end
    for id = poly_index
        M_poly = [id*ones(d,iLambda-1) indeces{1}];
        M = [M; M_poly];
    end
    indeces{iLambda} = M;                                                   % all index combinations of order iLambda
end

for iFile=1:K
disp(['Dataset_',num2str(iFile)])
%% Upload data
clear Input Output 
fileName = [num2str(iFile),dataset];
load(fileName);
if difFormat
    Input  = data_res(:,input_i)./normC;
    Output = data_res(:,output_i)./normC;
else
    Input  = fileData(:,input_i)./normC;
    Output = fileData(:,output_i)./normC;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the batch of input vectors
iNarx = 0;                                                                  % batch index of the input vector in AR model
timesNarx = [t_0:t_0+T];
for t=timesNarx
    iNarx = iNarx + 1;
    switch direction
        case 'Backward'
            y_narx(iNarx,1) =  Output(t);
            for it=1:lambda
                x_narx(it,iNarx) = delta_operator(it,Output,t,dT,'Forward');
            end
        case 'Forward'
            y_narx(iNarx,1)  = delta_operator(lambda,Output,t,dT,direction);
            x_narx(1,iNarx)  = Output(t);
            for it=2:lambda
                x_narx(it,iNarx) = delta_operator(it-1,Output,t,dT,direction);
            end
    end   
    x_narx(lambda+1,iNarx)  = Input(t-1);
    for it=lambda+2:2*lambda
        x_narx(it,iNarx)      = delta_operator(it-lambda-1,Input,t-1,dT,'Forward');
    end
 
end
nNarx = iNarx;                                                              % length of NARX input batch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the batch of delta-domain regressors
iTerm = 0;
if (iFile > 1)
    for iLambda = 1:lambda_p                                                  % iLambda - order of polynomial term
        for k = 1:size(indeces{iLambda},1)                                  % get index combination for the regressor in the sum
            iTerm = iTerm + 1;                                              % term index in the polynomial
            for iNarx = 1:nNarx                                             % going throught the full NARX batch
                term(iNarx,iTerm) = regressor(x_narx(:,iNarx),...
                                              indeces{iLambda}(k,:));       % compute the regressor (numeric)
            end
        end
    end
else
    for iLambda = 1:lambda_p                                                  % iLambda - order of polynomial term
        for k = 1:size(indeces{iLambda},1)                                  % get index combination for the regressor in the sum
            iTerm = iTerm + 1;                                              % term index in the polynomial
            for iNarx = 1:nNarx                                             % going throught the full NARX batch
                term(iNarx,iTerm) = regressor(x_narx(:,iNarx),...
                                              indeces{iLambda}(k,:));       % compute the regressor (numeric)
            end
            symb_term{iTerm} = strcat(x_str{indeces{iLambda}(k,:)});        % dictionary of regressors (symbolic)
        end
    end
end
iTerm = iTerm + 1;
term(:,iTerm) = 1;
symb_term{iTerm} = sym('c');
disp('Dictionary complete')
fileName = [folderName,'/dict_',dataset,num2str(iFile),'.mat'];
save(fileName, 'term','x_narx','y_narx','nNarx','t_0','-v7.3');
clear term x_narx y_narx
end
%% Save meta parameters
dictFolder = folderName;                                                    % folder from which I take dictionaties
nTerms = iTerm;                                                             % total number of regressors in the polynomial
dict_terms = [1:nTerms];                                                    % dictionary of all terms
fileMeta = ['Meta_delta_',dataset];
save(fileMeta,'dictFolder','nTerms','nNarx','symb_term','x_str','y_str','dict_terms','indeces','lambda','d','K','normC','-v7.3');   % save metadata
fileMeta = [folderName,'/Meta_delta_',dataset];
save(fileMeta,'dictFolder','nTerms','nNarx','symb_term','x_str','y_str','dict_terms','indeces','lambda','d','K','normC','-v7.3');
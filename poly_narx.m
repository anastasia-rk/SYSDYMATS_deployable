function poly_narx(foamset,dataset,normalise)
% Creates dictionaries of polynomial regressors from the time-series data
local_init;
%% Set up global parameters
difFormat = false;
switch foamset
    case 'foam_2010'
        input_i  = 3;                                                       % input column index
        output_i = 2;                                                       % output column index
        K        = 10;                                                      % number of datasets
        dT       = 0.02;
        normC    = 400;
    case 'foam_2019'
        input_i  = 2;                                                       % input column index
        output_i = 3;                                                       % output column index
        K        = 12;                                                      % number of datasets
        dT       = 2.3438e-04;
        normC    = 100;
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
        folder = 'dictionaries';                                             % specify category where to save files
        normC = 1;
        center = false;
    case 'Yes'
        folder = 'dictionaries_norm';                                        % specify category where to save files
        normC = 1;
        center = true;
end
% foamset = 'foam_2010'; % 'foam_2019'
% dataset = 'C'; %  'D'; %                                                    % name of dataset
addpath(['..\SYSDYMATS_data\',foamset])
iFile   = 1;                                                                % id of the sample
% Length of input and output lags
n_u     = 4;                                                                % input signal lag length
n_y     = 4;                                                                % output signal lag length
d       = n_y + n_u;                                                        % size of input vector x
lambda  = 3;                                                                % order of polynomial
% a       = sym('x_',[1 d]);                                                % associated symbolic vector
names = {'set','ny','nu'};                                                  % names used to define results folder name (no more than 3).
folderName = make_folder(folder,names,dataset,n_y,n_u);                     % create results folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create string array of input vector
iStr = 1;
if n_y~=0
   for it=n_y:-1:1
      x_str{iStr} = ['y(t-',num2str(it),')'];
      iStr = iStr + 1;
   end
end
for it=n_u:-1:1
      x_str{iStr} = ['u(t-',num2str(it),')'];
      iStr = iStr + 1;
end
y_str = 'y(t)';
% x_str{iStr} = ['u(t)'];   
%% Identify difference in lag
df = n_u - n_y;                                                           % difference between lags
switch sign(df)
    case 0
%         disp('zero')
        t_0 = n_u+1;
    case -1
%         disp('negative')
        t_0 = n_y+1;
    case 1
%         disp('positive')
        t_0 = n_u+1;
end
% t_0 = 500 + t_0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create sum index permutations
indeces{1} = [1:d]';
for iLambda=2:lambda                                                        % iLambda - order of polynomial term
    initial_m = permn(1:d,iLambda);                                         % get all permutations with repetition
    M = initial_m;                                                          % M - temporary set of permutations
    for j=iLambda:-1:2
        ind = find(M(:,j)>=M(:,j-1));                                       % sort out only increasing indeces
        M = M(ind,:);                                                       % update set
        clear ind
    end
    indeces{iLambda} = M;                                                   % all index combinations of order iLambda
end
for iFile=1:K
    if disFlag
       disp(['Dataset_',num2str(iFile)])
    end
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
if center
   meanIn   = mean(Input);
   Input  = Input - meanIn;
   meanOut   = mean(Output);
   Output  = Output - meanOut;
end
T   = min(10000,length(Input)-n_u);                                         %length(Input); % length of the observation sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the batch of input vectors
iNarx = 0;                                                                  % batch index of the input vector in AR model
timesNarx = [t_0:T];
for t=timesNarx
    iNarx = iNarx + 1;
    if n_y == 0
    x_narx(:,iNarx) = [Input(t-n_u:t,1)]; %                                 % NARX input
    else
        x_narx(:,iNarx) = [Output(t-n_y:t-1,1); Input(t-n_u:t-1,1)]; %      % NARX input
    end 
end
nNarx = iNarx;                                                              % length of NARX input batch
y_narx(:,:) = Output(timesNarx);                                            % NARX output
u_narx(:,:) = Input(timesNarx); % for calling the input from the inline function in MPO - no delay is u(t-1) called from the function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the batch of regressors and compute all polynomial terms
iTerm = 0;
if (iFile > 1)
    for iLambda = 1:lambda                                                  % iLambda - order of polynomial term
        for k = 1:size(indeces{iLambda},1)                                  % get index combination for the regressor in the sum
            iTerm = iTerm + 1;                                              % term index in the polynomial
            for iNarx = 1:nNarx                                             % going throught the full NARX batch
                term(iNarx,iTerm) = regressor(x_narx(:,iNarx),...
                                              indeces{iLambda}(k,:));       % compute the regressor (numeric)
            end
        end
    end
else
    for iLambda = 1:lambda                                                  % iLambda - order of polynomial term
        for k = 1:size(indeces{iLambda},1)                                  % get index combination for the regressor in the sum
            iTerm = iTerm + 1;                                              % term index in the polynomial
            for iNarx = 1:nNarx                                             % going throught the full NARX batch
                term(iNarx,iTerm) = regressor(x_narx(:,iNarx),...
                                              indeces{iLambda}(k,:));       % compute the regressor (numeric)
            end
            symb_term{iTerm} = strcat(x_str{indeces{iLambda}(k,:)});        % dictionary of regressors (symbolic)
            if size(indeces{iLambda},2) > 1
                inline_term{iTerm} = [sprintf('%s*',x_str{indeces{iLambda}(k,1:end-1)}),x_str{indeces{iLambda}(k,end)}];
            else
                inline_term{iTerm} = x_str{indeces{iLambda}(k,1)};
            end
        end
    end
end
iTerm = iTerm + 1;
term(:,iTerm) = 1;
symb_term{iTerm} = sym('c');
inline_term{iTerm} = '1';
if disFlag
    disp('Dictionary complete')
end
fileName = [folderName,'/dict_',dataset,num2str(iFile),'.mat'];
save(fileName, 'term','u_narx','x_narx','y_narx','nNarx','t_0','-v7.3');
clear term x_narx y_narx u_narx y u y_kf u_kf
end                                                                         % end loop over files
dictFolder = folderName;                                                    % folder from which I take dictionaties
nTerms = iTerm;                                                             % total number of regressors in the polynomial
dict_terms = [1:nTerms];                                                    % dictionary of all terms
fileMeta = ['Meta_',dataset];
fileMetaLocal = [folderName,'/Meta_',dataset];
save(fileMeta,'dictFolder','nTerms','nNarx','x_str','y_str','symb_term','inline_term','dict_terms','indeces','n_y','n_u','K','normC','center','-v7.3');
save(fileMetaLocal,'dictFolder','nTerms','nNarx','x_str','y_str','symb_term','inline_term','dict_terms','indeces','n_y','n_u','K','normC','center','-v7.3');   
        
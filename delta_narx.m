function delta_narx(foamset,dataset,direction,normalise)
% Creates dictionaries of delta-domain regressors from the time-series data
local_init;
visFlag = true;
%% Set up global parameters
difFormat = false;
switch foamset
    case 'foam_2010'
        input_i  = 3;                                                       % input column index
        output_i = 2;                                                       % output column index
        K        = 10;                                                      % number of datasets
        normC    = 400;
        dT       = 0.02;
        % for KF derivative estimation - process noise variances
        sig2_y = 1000;
        sig2_u = 1000;
    case 'foam_2019'
        input_i  = 2;                                                       % input column index
        output_i = 3;                                                       % output column index
        K        = 12;                                                      % number of datasets
        normC    = 100;
        dT       = 2.3438e-04;
        difFormat = true;
        % for KF derivative estimation - process noise variances
        sig2_y = 10^10;
        sig2_u = 10^10;
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
        % for KF derivative estimation - process noise variances
        sig2_y = 10^8;
        sig2_u = 10^11;
        difFormat = true;
    case 'sweep_2020'
        input_i  = 2;                                                       % input column index
        output_i = 3;                                                       % output column index
        dT       = 9/51200;                                                 % sampling time
        K        = 6;
        normC    = 1;
        % for KF derivative estimation - process noise variances
        sig2_y = 10^10;
        sig2_u = 10^12;
        normC = 1;
end
folder = 'delta_dictionaries';                                              % specify category where to save files
switch direction
    case 'Backward'
        folder = [folder,'_b'];
    case 'Forward'
        folder = [folder,'_f'];
end
switch normalise
    case 'No'
        normC = 1;
    case 'Yes'
        folder = [folder,'_norm'];                                          % normalised data dictionaries
end

addpath(['../SYSDYMATS_data/',foamset])
% Length of input and output lags
lambda  = 2;                                                                % order of delta operator
lambda_p = 3;    
n_u = 1+lambda;
n_y = lambda;% order of polynomial
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
                        for it=lambda+2:2*lambda+1
                            x_str{it} = ['\delta^',num2str(it-lambda-1),' u(t)'];
                        end
                case 'Forward'
                        folder = [folder,'_f'];
                        x_str{1} = 'y(t)';
                        for it=2:lambda
                            x_str{it} = ['\delta^',num2str(it-1),' y(t)'];
                        end
                        x_str{lambda+1} = 'u(t)';
                        for it=lambda+2:2*lambda+1
                            x_str{it} = ['\delta^',num2str(it-lambda-1),' u(t)'];
                        end
                        poly_index = [1 lambda+1];
                        y_str = ['\delta^',num2str(lambda),' y(t)'];
                        
end
        d = it;
t_0 = lambda+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Smoothers for discrete derivative estimation
plotInt = 1:500;
sgolayOrder = 15;
sgolayFrame = min(25,sgolayOrder*2+1);
% Savitzky-Golay - polynomial moving window smoother that preserves signal
% shape with high order polynomials
[b,g] = sgolay(sgolayOrder,sgolayFrame);
% KF - measured input and output cnsidered as noiseless measurements
A = eye(lambda+1);
for i=1:lambda+1
    for j=i+1:lambda+1
        A(i,j) = dT^(j-i)/factorial(j-i);
    end
end
% A = [zeros(lambda,1) triu(ones(lambda+1))*dT; zeros(1,lambda+1)];
B = zeros(lambda+1,1); u = 0; % no deterministic input
G = [zeros(lambda,1); 1];
C = [1 zeros(1,lambda)];
Q_y = sig2_y*G*G';
Q_u = sig2_u*G*G';
R = 0.0001;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main loop
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T   = min(10000,length(Input)-n_u);                                         % length of the observation sequence
iNarx = 0;                                                                  % batch index of the input vector in AR model
timesNarx = [t_0:T];
y = Output(timesNarx);
u = Input(timesNarx-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Direct delta-op computation
 for t=lambda+1:length(timesNarx)
        for iLambda=1:lambda
            delta_y_b(iLambda,t) = delta_operator(iLambda,y,t,dT,'Backward'); 
            delta_u_b(iLambda,t) = delta_operator(iLambda,u,t,dT,'Backward'); 
        end
 end
  for t=1:length(timesNarx)-lambda
        for iLambda=1:lambda
            delta_y_f(iLambda,t) = delta_operator(iLambda,y,t,dT,'Forward'); 
            delta_u_f(iLambda,t) = delta_operator(iLambda,u,t,dT,'Forward'); 
        end
 end
%% Savitzky-Golay filter for smoothing derivatives
        for it = 0:lambda
            dy(:,it+1) = conv(y, factorial(it)/(-dT)^it * g(:,it+1), 'same');
        end
        
        for it = 0:lambda
            du(:,it+1) = conv(u, factorial(it)/(-dT)^it * g(:,it+1), 'same');
        end
%% Kalman filter dor derivative estimation
x_y = zeros(lambda+1,1);
P_y = eye(lambda+1);
x_u = zeros(lambda+1,1);
P_u = eye(lambda+1);
for t=1:length(timesNarx)
     [x_y,P_y,~] = kf(y(t),x_y,0,P_y,A,B,C,Q_y,R);
     [x_u,P_u,~] = kf(u(t),x_u,0,P_u,A,B,C,Q_u,R);
     y_kf(:,t) = x_y;
     P_y_kf{t} = P_y;
     u_kf(:,t) = x_u;
     P_u_kf{t} = P_u;
end
for t=length(timesNarx):-1:1
     [x_y,P_y] = rts(x_y,y_kf(:,t),P_y,P_y_kf{t},0,A,B,Q_y);
     [x_u,P_u] = rts(x_u,u_kf(:,t),P_u,P_u_kf{t},0,A,B,Q_u);
     y_ks(:,t) = x_y;
     P_y_ks{t} = P_y;
     u_ks(:,t) = x_u;
     P_u_ks{t} = P_u;
end
%% Plots
figName = [dataset,num2str(iFile)];
L2 = lambda + 1;
fig(figName,visFlag); 
subplot(2,L2,1);
plot(y(plotInt)); hold on;
plot(dy(plotInt,1)); hold on;
plot(y_kf(1,plotInt)); hold on;
plot(y_ks(1,plotInt)); hold on;
xlabel('Index');ylabel('$y(t)$');
legend('Measured','SGF','KF','RTS');
for iLambda = 1:lambda
    subplot(2,L2,iLambda+1);
    plot(lambda+plotInt,delta_y_b(iLambda,plotInt)); hold on;
    plot(delta_y_f(iLambda,plotInt)); hold on;
    plot(dy(plotInt,iLambda+1)); hold on;
    plot(y_kf(iLambda+1,plotInt)); hold on;
    plot(y_ks(iLambda+1,plotInt)); hold on;
    xlabel('Index');ylabel(['$\delta^{',num2str(iLambda),'}y(t)$']);
%     legend('$\delta$ forward','SGF','KF','RTS');
    legend('$\delta$ backward','$\delta$ forward','SGF','KF','RTS');
end
subplot(2,L2,lambda+2);
plot(u(plotInt)); hold on;
plot(du(plotInt,1)); hold on;
plot(u_kf(1,plotInt)); hold on;
plot(u_ks(1,plotInt)); hold on;
xlabel('Index');ylabel('$u(t)$');
legend('Measured','SGF','KF','RTS');
for iLambda = 1:lambda
    subplot(2,L2,iLambda+lambda+2);
    plot(lambda+plotInt,delta_u_b(iLambda,plotInt)); hold on;
    plot(delta_u_f(iLambda,plotInt)); hold on;
    plot(du(plotInt,iLambda+1)); hold on;
    plot(u_kf(iLambda+1,plotInt)); hold on;
    plot(u_ks(iLambda+1,plotInt)); hold on;
    xlabel('Index');ylabel(['$\delta^{',num2str(iLambda),'}u(t)$']);
%     legend('$\delta$ forward','SGF','KF','RTS');
    legend('$\delta$ backward','$\delta$ forward','SGF','KF','RTS');
end
tikzName = [folderName,'/','derivatives_',dataset,num2str(iFile),'.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '15cm', 'width','15cm','checkForUpdates',false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the batch of input vectors
switch direction
       case 'Backward'
            y_narx = y_ks(1,:)';
            x_narx = [y_ks(2:lambda+1,:); u_ks];
       case 'Forward'
            y_narx = y_ks(end,:)';
            x_narx = [y_ks(1:lambda,:); u_ks];
end
nNarx = length(y_narx);                                                     % length of NARX input batch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the batch of delta-domain regressors
iTerm = 0;
if (iFile > 1)
    for iLambda = 1:lambda_p                                                % iLambda - order of polynomial term
        for k = 1:size(indeces{iLambda},1)                                  % get index combination for the regressor in the sum
            iTerm = iTerm + 1;                                              % term index in the polynomial
            for iNarx = 1:nNarx                                             % going throught the full NARX batch
                term(iNarx,iTerm) = regressor(x_narx(:,iNarx),...
                                              indeces{iLambda}(k,:));       % compute the regressor (numeric)
            end
        end
    end
else
    for iLambda = 1:lambda_p                                                % iLambda - order of polynomial term
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
if disFlag
    disp('Dictionary complete')
end
fileName = [folderName,'/dict_',dataset,num2str(iFile),'.mat'];
save(fileName, 'term','x_narx','y_narx','nNarx','t_0','-v7.3');
clear term x_narx y_narx
end
%% Save meta parameters
dictFolder = folderName;                                                    % folder from which I take dictionaties
nTerms = iTerm;                                                             % total number of regressors in the polynomial
dict_terms = [1:nTerms];                                                    % dictionary of all terms
switch direction
       case 'Backward'
            fileMeta = ['Meta_delta_',dataset,'_b'];
       case 'Forward'
            fileMeta = ['Meta_delta_',dataset,'_f'];
end
save(fileMeta,'dictFolder','nTerms','nNarx','symb_term','x_str','y_str','dict_terms','indeces','lambda','d','K','normC','-v7.3');   % save metadata
fileMeta = [folderName,'/Meta_delta_',dataset];
save(fileMeta,'dictFolder','nTerms','nNarx','symb_term','x_str','y_str','dict_terms','indeces','lambda','d','K','normC','-v7.3');
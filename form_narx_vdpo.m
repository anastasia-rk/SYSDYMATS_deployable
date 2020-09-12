local_init;
% Creates dictionaries of delta-domain regressors from the time-series data
%% Set up global parameters
visFlag = true;
foamset     ='vdpo';
input_i     = 2;                                                       % input column index
output_i = 3;                                                       % output column index
K        = 10;                                                      % number of datasets
normC    = 1;
dataset = 'V';
addpath('..\SYSDYMATS_data\',foamset)
load('External_parameters_V');
folder = 'dictionaries';
normalise = questdlg('Scale the data?', ...
        'Normalisation');
switch normalise
    case 'No'                                     
        normC = 1;
    case 'Yes'
        folder = [folder,'_norm'];
        normC = 4;
end
regressors = questdlg('Select the domain of regressors', ...
    'Domain choice',...
	'Shift','Delta','');
switch regressors
    case 'Shift'
        n_u     = 4;                                                        % input signal lag length
        n_y     = 4;                                                        % output signal lag length
        d       = n_y + n_u;                                                % size of input vector x
        lambda_p  = 3;                                                      % order of polynomial
        names = {'set','ny','nu'};                                          % names used to define results folder name (no more than 3).
        folder = ['../SYSDYMATS_dictionaries/',folder];
        folderName = make_folder(folder,names,dataset,n_y,n_u);             % create results folder
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

    case 'Delta'
        folder = ['delta_',folder];
        direction = questdlg('Type of delta operator', ...
                    'Causality',...
                    'Forward','Backward','');
        lambda  = 2;                                                        % order of delta operator
        lambda_p = 3;                                                       % order of polynomial
        n_u = 1+lambda;
        n_y = lambda;
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
        names = {'set','lambda'}; 
        folder = ['../SYSDYMATS_dictionaries/',folder];
        folderName = make_folder(folder,names,dataset,lambda);             % create results folder

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Smoothers for discrete derivative estimation
dT  = 1/sampFr;
plotInt = 1:500;
if exist('lambda') == 1
    sgolayOrder = 9;
    sgolayFrame = min(45,sgolayOrder*2+1);
    % Savitzky-Golay - polynomial moving window smoother that preserves signal
    % shape with high order polynomials
    [b,g] = sgolay(sgolayOrder,sgolayFrame);
    % KF - measured input and output cnsidered as noiseless measurements
    sig2_y = 1;
    sig2_u = 4;
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create index permutations for sums
% Intial combination accounts for all linear regressors
indeces{1} = [1:d]';
switch regressors
    case 'Shift'
    % With self-products
    for iLambda=2:lambda_p                                                      % iLambda - order of polynomial term
    initial_m = permn(1:d,iLambda);                                         % get all permutations with repetition
    M = initial_m;                                                          % M - temporary set of permutations
    for j=iLambda:-1:2
        ind = find(M(:,j)>=M(:,j-1));                                       % sort out only increasing indeces
        M = M(ind,:);                                                       % update set
        clear ind
    end
    indeces{iLambda} = M;                                                   % all index combinations of order iLambda
end
    case 'Delta'
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
end 
%% Main loop
for iFile=1:K
    if disFlag
       disp(['Dataset_',num2str(iFile)])
    end
% Upload data
clear Input Output 
fileName = [num2str(iFile),dataset];
load(fileName);
Input  = fileData(:,input_i)./normC;
Output = fileData(:,output_i)./normC;
t_0 = 100;
T   = min(10000,length(Input)-n_u); %length(Input); % length of the observation sequence
iNarx = 0;                                                                  % batch index of the input vector in AR model
timesNarx = [t_0:T];
y = Output(timesNarx);
u = Input(timesNarx-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Direct delta-op computation
if exist('lambda') == 1
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
    legend('$\delta$ backward','$\delta$ forward','SGF','KF','RTS');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the batch of input vectors
switch regressors
    case 'Shift'
        for t=timesNarx
            iNarx = iNarx + 1;
            if n_y == 0
                x_narx(:,iNarx) = [Input(t-n_u:t,1)]; %                     % NARX input
            else
                x_narx(:,iNarx) = [Output(t-n_y:t-1,1);Input(t-n_u:t-1,1)]; % NARX input
            end 
        end
        y_narx(:,:) = Output(timesNarx);                                    % NARX output
    case 'Delta'
            switch direction
                case 'Backward'
                    y_narx = y_ks(1,:)';
                    x_narx = [y_ks(2:lambda+1,:); u_ks];
                case 'Forward'
                    y_narx = y_ks(end,:)';
                    x_narx = [y_ks(1:lambda,:); u_ks];
              end   
end
nNarx = length(y_narx);  
u_narx(:,:) = Input(timesNarx); % for calling the input from the inline function in MPO - no delay is u(t-1) called from the function
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
clear term x_narx y_narx u_narx y u y_kf u_kfend
end
%%
% end loop over files
dictFolder = folderName;                                                    % folder from which I take dictionaties
nTerms = iTerm;                                                             % total number of regressors in the polynomial
dict_terms = [1:nTerms];                                                    % dictionary of all terms
% save metadata
switch regressors
    case 'Shift'
        fileMeta = ['Meta_',dataset];
        fileMetaLocal = [folderName,'/Meta_',dataset];
        save(fileMeta,'dictFolder','nTerms','nNarx','x_str','y_str','symb_term','inline_term','dict_terms','indeces','lambda_p','n_y','n_u','K','normC','-v7.3');
        save(fileMetaLocal,'dictFolder','nTerms','nNarx','x_str','y_str','symb_term','inline_term','dict_terms','indeces','lambda_p','n_y','n_u','K','normC','-v7.3');   
    case 'Delta'
        fileMeta = ['Meta_delta_',dataset];
        fileMetaLocal = [folderName,'/Meta_delta_',dataset];
        save(fileMeta,'dictFolder','nTerms','nNarx','symb_term','inline_term','x_str','y_str','dict_terms','indeces','lambda','lambda_p','d','K','normC','-v7.3');   
        save(fileMetaLocal,'dictFolder','nTerms','nNarx','symb_term','inline_term','x_str','y_str','dict_terms','indeces','lambda_p','d','K','normC','-v7.3');
end

%% Add mylibs to and subdirectories to path in local directory
% Requires that mylibs is located in the parent folder.
%userpath reset; startup;
LibFolder = 'mylibs';                                                       % name of library name
pathCell = regexp(path, pathsep, 'split')';                                 % array of directories in Matlab path
index_path = find(cell2mat(strfind(pathCell,LibFolder)));                   % find whether path contains mylibs
LibFolder = ['../',LibFolder];
if isempty(index_path)
    addpath(LibFolder)
    newPath = true;
else
    myPath = pathCell(index_path);                                          % see all subridectories in the path
    newPath = false;
end
%% Go through the subdirectories
f = dir(LibFolder);
git = 4;                                                                    % git offset
N = length(f);
for iDir=git:N
    if f(iDir).isdir
       if newPath || ~any(cell2mat(strfind(myPath,f(iDir).name)))
        addpath([LibFolder,'/',f(iDir).name]);
       end
    end
end
%% Figure and interpreter setup
clear all; my_init;
visFlag = 'On'; %% or 'Off'
disFlag = true; % set to false to prevent from displaying updates in work space
%% CVX setup for SHARC - have to call every time
% addpath('../MATLAB/cvx');
% cvx_setup  cvx_license.dat
%% Optimisation and parfor settings
%% before this function open matlab and run configCluster;
function submit_matlab_parallel

cd path_working_directory;
c = parcluster;
c.AdditionalProperties.EmailAddress = 'a.kadochnikova@sheffield.ac.uk';
% Configure runtime e.g. 40 minutes
c.AdditionalProperties.WallTime = '60:00:00';
% Configure rmem per process e.g. 4 Gb
c.AdditionalProperties.AdditionalSubmitArgs = ' -l rmem=4G';
% Parallel_example.m contains the parfor loop, no_of_cores < 31
j = c.batch(@parallel_example, 1, {}, 'Pool', 23);
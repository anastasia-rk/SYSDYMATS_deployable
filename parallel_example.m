function time = parallel_example
cd path_working_directory;
outfile = ['output.txt'];
fileID = fopen(outfile, 'w');
%disp('Parallel time')
tic
n = 200;
A = 500;
a = zeros(n);
parfor i = 1:n
    a(i) = max(abs(eig(rand(A))));
end
time = toc;
fprintf(fileID, '%d', time);
fclose(fileID);
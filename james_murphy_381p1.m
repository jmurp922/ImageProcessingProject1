
% HW 1  
%-----------------------------PART 1----------------------------------%
mat_a = randi([0, 10], [100, 100]);
mat_b = randi([0, 10], [100, 100]);

%Declare the table for the graph data
graph_data(1, 100) = 100;

% find the size of matrix as columns & rows
% m is number of rows
% n is number of columns

[m, n] = size(mat_a);

%Zero out & initalize a matrix of size m x n which is equal to the size of
%matrix a

mat_c = zeros(m, n);

%tic is the start of the stopwatch timer in matlab
tic;

for i = 1 : m
    for j = 1 : n
   
        %to find the i, j entry of the euclidian distance you take the row 
        %of matrix a from i to m and subtract it
        %from the row of matrix b from j to m and then square it, sum it and then
        %take the square root of it according to the formula given
        mat_c(i, j) = sqrt(sum((mat_a(i, :) - mat_b(j, :)).^2));
        
        %increment j
        j = j + 1;
    end
    %increment i
    i = i + 1;
    % Toc stops the timer, store it in variable end_time so we can update
    % the graph on the elapsed time.
    end_time_pt1 = toc;
    %-----------------------------PART 2----------------------------------%
   
    tic;
    
    % Initalize a variable matrix d and subtract a and b and store it in d.
    mat_d = mat_a - mat_b;

    % take the sum of the matrix 
    euc_dist = sum(mat_d);

    % square the matrix 
    euc_dist = (euc_dist.^2);

    % take the square root of the matrix to get the euclidian distance
    euc_dist = sqrt(euc_dist);

    end_time_pt2 = toc;

end

    
v0 = graph_data(1,:);
v1 = graph_data(2,:);
v2 = graph_data(3,:);

%Creates a window on the screen
figure

%Plot part 1
plot(v0,v1);

%Don't delete existing plots
hold on

%Plot part 2
plot(v0,v2);

%Name labels for x and y axis, and the legend to designate which operation
%is being done.
legend('Loop','Matrix Operations')
xlabel("Columns")
ylabel('Running Time (Seconds)')


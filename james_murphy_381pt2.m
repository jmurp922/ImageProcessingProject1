%HW 1 PART 2

% Create an arbitrary max integer for the loop to run through. 
max_num = 100;  

% Create an arbitrary integer to start at, in this case we will be doing 1 to max_num
% (100) iterations

start_num = 1;

%iteration & operation solution, loops from 1 -> 100

%matrix operation solution is at the end of hte loop.
for start_num = 1 : max_num
    %-----------------------------PART 1----------------------------------%
    %tic is the start of the stopwatch timer in matlab
    %tic is placed inside the loop so that we can keep track of the time
    %for each iteration up to max_num
    tic
    
    % Create a new variable graphdata. graphData 
    % goes from 1 -> 100 which is the maximum number we allowed
    
    graph_data(1, start_num) = start_num;
    
    %disp("Graph data is : " + graphData);
    
    %Create a matrix with a random number with dimensions of the number of  
    %iterations currently being ran through. (ex : on iteration 1 -> 1, 
    % 2 -> 2, .... , 100 -> 100 (Max num).  
    
    rand_matrix = randi([0,10],[start_num,start_num]);
    
    % Need 2 matrices of the same dimension to multiply to solve the
    % problem, so mat_a and mat_b are zeroed out and are made to be the
    % same dimension..
    
    mat_a = zeros(start_num,start_num);
    
    mat_b = zeros(start_num,start_num);
    
    %Create the average of the graph table and zero it out from 1 -> the
    %current loops iteration 
    avg_table = zeros(1, start_num);
    
    for index = 1 : start_num
        %Initialize the result to 0
        result = 0;
        
        % Another for loop to calculate the result. Take the result and add
        % it to the random matrix starting from the inner_loop value (rows)
        % and the total number of columns allowed in our execution, start_num 
        %(1 -> 100)
         for inner_loop = 1 : start_num
           result = result + rand_matrix(inner_loop, start_num);
           
           %increment the loop
           
           inner_loop = inner_loop + 1;
         end
         
         % Take the value of result after the inner loop and divide it by
         % the start number and then put that result in the average table
         % at the current index in the loop. 
         
         result = result / start_num;
         avg_table(index) = result;
         
         %increment the loop
         
         index = index + 1;
    end
    
    % Now another loop, we have to get the covariance
    
    % In this loop we're getting the Covariance matrix. Since this is the
    % nested loop solution we start loop i and j as m & n. M is the number of
    % rows, n is the number of columns. We need to iterate accross the
    % whole matrix and set the value at i, j in matrix a to the result we 
    % derive from our covariance formula
    for m = 1 : start_num
        for n =1 : start_num
            % This loop 3rd nested loop is used to calculate the covariance
            % Declare the result that we operate in on the covariance loop
             result = 0;
            for covariance_loop = 1 : start_num
               %Take the result and add it with the current index of the
               %covariance_loop and the outmost row loop (m). Then,
               %subtract that from the avg_table (mean) row and multiply it
               %by the same equation we just did but swapping the row with columns 
               %(swap m and n)
               
               result = result + (rand_matrix(covariance_loop, m) - avg_table(m) * (rand_matrix(covariance_loop, n) - avg_table(n)));
               
               %Then, store that result into the i, j pair of our current
               %matrix.
               mat_a(m, n) = result;
               
               %increment the loop
               covariance_loop = covariance_loop + 1;
               
            end 
            
            %increment the loop
            n = n + 1;
        end
        
        %increment the loop
        m = m +1;
    end
    
    %Now we need to calculate the standard deviation
    
    %Take the diagonal of our matrix
    standard_deviation = diag(mat_a);
    
    %Tranpose the diagonal of the matrix
    standard_deviation = transpose(standard_deviation);
    
    %Take the square root of the transpose of the diagonal 
    standard_deviation = sqrt(standard_deviation);
    
    %M = rows, N = cols
    
    for m = 1:start_num
        for n = 1:start_num
            
            %We take the zerod out matrix b, and take the result of the i,
            %j pair of covariance matrix i and divide it by the standard
            %deviation using the rows variable m

            mat_b(m,n) = mat_a(m, n) / standard_deviation(m);
            
            % After that we take the value we jsut stored at the i, j pair
            %and divide it by the standard deviation.
            
            mat_b(m,n) = mat_b(m, n) / standard_deviation(n);
            
            %increment the loop
            n = n + 1;
        end
        
        %increment the loop
        m = m + 1;
        
    end
    
    % mat_b is now the correlation matrix
    
    % Toc stops the timer, store it in variable end_time so we can update
    % the graph on the elapsed time.
    
    timep1 = toc; 
    
    %See individual time iterations 
    %disp("End Time at iteration : " + start_num + " is : " + end_time);
    
    %-----------------------------PART 2----------------------------------%
    tic 

    % Already calculated the mean matrix so we can cheat and use it.
    avg_mat = avg_table;
    
    % take the random matrix and subtract it from the mean we calculated
    current_mat = rand_matrix - avg_mat
    
    % Now we need the covariance matrix
    % Take the transpose of the matrix we just calculated
    trans_mat = transpose(current_mat);
    
    % Now take the transpose matrix and * by the current_matrix and divide
    % it by the starting number - 1.
    
    covariance_matrix = trans_mat * (current_mat / (start_num - 1));
    
    % Now we need to assemble the D x D matrix as asked for in the main
    % part.
    
    %We need the diagonal of the covariance matrix.
    D = diag(covariance_matrix);
    
    %Take the square root of itself.
    D = sqrt(D);
    
    % Multiply the matrix by the transpose
    D = D * transpose(D);
    
    % Take the covariancematrix and divide it by D
    % (./ means Element-wise right -> division)
    
    correlation_matrix = covariance_matrix ./ D;    
    
    timep2 = toc; 
    
    % Part 1 is stored here each loop iteration
    graph_data(2,start_num) = timep1;
    
    % Part 2 is stored here each loop iteration
    graph_data(3,start_num) = timep2;
    
end

%disp(avg_table);

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

    
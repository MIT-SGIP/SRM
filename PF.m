%% Description

%{

    Identify which input points lie on the pareto front based on the 
    objective of minimizing max acceleration and max length.

    Input:
    - inputs: a nx2 matrix with design point data

    Output:
    - pareto: a mx2 matrix (m <= n) with points from "inputs" that lie on
    the pareto front that minimizes the first column and minimizes the
    second column

%}

%% Function

function pareto = PF(inputs)

    % if dominated = 1 for a given entry, it should be removed
    dominated = zeros(length(inputs), 1);

    for i = 1:length(inputs)

        for j = i+1:length(inputs)

            % check whether one element dominates the other
            if inputs(i, 1) > inputs(j, 1) && inputs(i, 2) > inputs(j, 2) 
                
                dominated(i) = 1;

            elseif inputs(i, 1) < inputs(j, 1) && inputs(i, 2) < inputs(j, 2)

                dominated(j) = 1;

            end

        end
        
    end

    % only include the non-dominated results 
    pareto = inputs(~dominated, :);

end
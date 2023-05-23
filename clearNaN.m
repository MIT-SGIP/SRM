%% Description

%{

    This function combs through the impact mechanics results array for one 
    impact case to remove rows of all NaNs. Specifically, this function
    identifies NaNs in the second column. 

    Inputs:
    - dataArr: array of all impact mechanics data, which has some rows of
    all NaN data.

    Outputs:
    - dataArrCorr: array of impact mechanics results without rows of NaNs.
    This data has the same number of columns as dataArr. dataArrCorr will
    have fewer rows than dataArr if NaNs are detected.

%}

%% Function

function dataArrCorr = clearNaN(dataArr)

    logicArr = isnan(dataArr(:,2));

    inds1 = (1:length(logicArr)).';

    % the indices of only NaN entries
    inds = find(inds1 .* logicArr);

    dataArrCorr = dataArr;

    % clear the NaN entries
    dataArrCorr(inds, :) = [];

end
%% Description

%{

    Preps the data from the team's Ansys impact simulation. The impact
    simulation outputs text files with 12 columns. Only columns 3, 5, 6, 8,
    9, 11, and 12 are used. See text files in "Input Data" folder for
    names and units.

    This script outputs a datafile called "PreppedData.mat" which when
    loaded in the structural response model, loads a "data" variable.

%}

clc;clear;close all;

warning off;

debug = 0;

%% Constants

addpath('Input Data');

filenames = ["dp0_66kg_0z_42.5y.txt", "dp5_66kg_2.5z_42.5y.txt", ...
    "dp10_66kg_5z_42.5y.txt", "dp16_66kg_7.5z_42.5y.txt", ...
    "dp21_66kg_10z_42.5y.txt"];

g = 9.81;

% temporary end times [s] 
endTimes = .25 * ones(1, length(filenames)); 

% initial velocities set by impact case [m/s]
initialVsZ = -2.5 * (0:1:4);
initialVsY = -42.5 * ones(1,length(filenames));

% data sample rate [s]
dtDes = 1e-4;

% initialize arrays
dtArr = [];
deltaVsZ = [];
deltaVsY = [];

%% Prep Data from each file

%{
 
    Data will contain the following format:
    1: t [s]
    2: accY [m/s2]
    3: PosYTip [m]
    4: PosZTip [m]
    5: PosYTail [m]
    6: accZ [m/s2]
    7: PosZTail [m]
    8: angles [deg]

%}

data = {};

for i=1:length(filenames)
    
    data1 = readtable(filenames(i));
    
    % convert tables to matrix
    dataArr = table2array(data1(:,[3,5,6,8,9,11,12]));
    
    % remove any NaNs from the results array
    dataArrCorr = clearNaN(dataArr);

    % find the angles for this 
    angles = atand((dataArrCorr(:,7) - dataArrCorr(:,4)) ./ (dataArrCorr(:,5) - dataArrCorr(:,3)));

    % append the angles to the end of this dataset
    dataAddition = [dataArrCorr, angles];

    % the current data's time vector
    tt = dataAddition(:,1);

    % the desired time vector (with more datapoints)
    ttDes = (0:dtDes:tt(end))';

    % spline interpolation
    dataAdditionSpline = zeros(length(ttDes), size(dataAddition,2));

    dataAdditionSpline(:,1) = ttDes;

    for j=2:size(dataAddition, 2)

        dataAdditionSpline(:,j) = spline(tt, dataAddition(:,j), ttDes);

    end

    [~, endIndTemp] = min(dataAdditionSpline(1:(end*.9),3));
    [~, endIndTemp2] = min(dataAddition(1:(end*.9),3));

    % figure 3-5 in Mike's thesis
    if i==1
        fig49 = figure;
        fig49.Position = [50, 50, 1000, 750];

        axialAccel222 = cosd(dataAddition(1:endIndTemp2,end)) .* dataAddition(1:endIndTemp2,2) + ...
            sind(dataAddition(1:endIndTemp2,end)) .* dataAddition(1:endIndTemp2,6);

        axialAccelSpline222 = cosd(dataAdditionSpline(1:endIndTemp,end)) .* dataAdditionSpline(1:endIndTemp,2) + ...
            sind(dataAdditionSpline(1:endIndTemp,end)) .* dataAdditionSpline(1:endIndTemp,6);
        
        hold on
        plot(ttDes(1:endIndTemp), axialAccelSpline222/g, 'LineWidth',2)
        scatter(tt(1:endIndTemp2), axialAccel222/g,400,'.k')
        hold off

        ax = gca;

        fs = 18;
        ax.FontSize = fs+2;
                xlabel('Time [s]', 'FontSize', fs+2)
        ylabel('Axial Acceleration [g]', 'FontSize', fs+2)

        [~, obj] = legend('Spline Interpolation', 'Ansys Output Data', ...
            'Location','best','FontSize', fs);

        objl = findobj(obj, 'type','patch');
        set(objl, 'Markersize', 25)
        
    end

    % add the current corrected data array to data storage
    data{i} = dataAdditionSpline;

    % estimated change in velocities [m/s]
    dvvy = cumtrapz(dataAddition(:,1), dataAddition(:,2));
    dvvz = cumtrapz(dataAddition(:,1), dataAddition(:,6));

    % find where the penetrator stops in the Y direction
    endInd = find((dvvy + initialVsY) > 0, 1);

    % the body doesn't rebound or stop in all sims. Add a check such that
    % if the body is still moving, use the sim's final time as the end time
    if length(endInd) > 0
        endTimes(i) =  dataAddition(endInd,1);
        % only get the change in velocity up to the cutoff time (when delta v Y is
        % 42.5 m/s)
        deltaVsY(1,i) = dvvy(endInd);
        deltaVsZ(1,i) = dvvz(endInd);
    else
        endTimes(i) = dataAddition(end,1);
        deltaVsY(1,i) = dvvy(end);
        deltaVsZ(1,i) = dvvz(end);
    end
    
    diffArrMax = max(diff(dataArrCorr(:,1)));
    dtArr = [dtArr,diffArrMax];

end

save("PreppedData.mat", 'data')

%% data for table 3-1 of Mike's thesis

% initialize data storage for table
datStore = zeros(5,4);

for i = 1:length(filenames)
    
    % get the load case of interest
    datt = data{i};

    % find the index of the maxmimum penetration distance
    [~, endInd] = min(datt(1:(end*.9),3));

    % find the end time 
    endTime = datt(endInd, 1);

    % tipping angle
    theta = datt(1:endInd, end);

    % Y and Z accelerations
    accZZ = datt(1:endInd, 6);
    accYY = datt(1:endInd, 2);

    % get the axial and lateral accelerations over this time
    lateralAccel = accZZ .* cosd(theta) - accYY .* sind(theta); 
    axialAccel = accYY .* cosd(theta) + accZZ .* sind(theta);

    % plot figures if debugging
    if debug
        figure
        hold on
        plot(datt(1:endInd,1), datt(1:endInd,3))
        hold off
        xlabel('t [s]')
        ylabel('penetration distance')
    
        figure
        hold on
        plot(datt(:,1), datt(:,3))
        hold off
        xlabel('t [s]')
        ylabel('penetration distance')
    
        figure
        hold on
        plot(datt(1:endInd,1), axialAccel/g)
        plot(datt(1:endInd,1), lateralAccel/g)
        hold off
        xlabel('t [s]')
        ylabel('acceleration [g]')
        legend('axial','lateral')
    end

    % get the maximum magnitude of axial or lateral acceleration

    % axial case
    maxMagAx = max(axialAccel);

    if abs(min(axialAccel)) > maxMagAx
        maxMagAx = min(axialAccel);  
    end

    avgAccelAx = mean(axialAccel);

    % lateral case
    maxMagLat = max(lateralAccel);

    if abs(min(lateralAccel)) > maxMagLat
        maxMagLat = min(lateralAccel);  
    end
    
    % store data for this load case
    datStore(i, :) = [endTime, avgAccelAx/g, maxMagAx/g, maxMagLat/g];

end

%% Save end times

endTimes = datStore(:,1);

save("ImpactDurations.mat", "endTimes");
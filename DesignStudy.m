%% Description

%{

    This script implements the design study evaluating different shock
    isolation design points. Each design point is run through all 5 impact
    cases described in Mike's thesis.

    dataPrep.m must be run or dataPrep.mat and ImpactDurations.mat must be 
    present to run this script.

%}

clc;clear;close all;

%% Constants 

% load the prepped data
load("PreppedData.mat");
load("ImpactDurations.mat")

% gravity constant [m/s2]
g = 9.81;  

% masses of objects [kg]
mBody = 66.60254149;
mFlare = 25.59132169;
mSeismometer = 2;   

% store masses for DP function
masses = [mBody, mFlare, mSeismometer];

% flags to show data for debugging
debug = 0;          % set to 1 to see each design undergoing each impact case
filterDebug = 0;    % set to 1 to see how the 

% initial velocities [m/s]
initialVsY = -42.5 * ones(1, 5);
initialVsZ = -2.5 * (0:1:4);

% threshold for when to stop simulation based on seismometer energy [J]
energyThreshold = 1e-2; 

%% Assemble the control parameters

rigidBodyStore = zeros(1, 3);

matProps = {};
model = @controlModel;

% design labelled with 1 is a crushable material
rigidInd = 0;

rigidResultStore(1,1:2) = DP(data, model, matProps, debug, ...
    filterDebug, rigidInd, masses, g, endTimes, energyThreshold, ...
    initialVsY, initialVsZ);

rigidResultStore(:,3) = ones(size(rigidResultStore, 1), 1) * rigidInd;

% resultStore is in form: max displacement [m], max axial acceleration
% [m/s2], index of material type

%% Assemble different paremeters to track between crushable models

disp("Starting crushable material models")

s = tic; 

% crushable plateau force [N]
plateauForce = 800:40:2000;

% space efficiency of shock isolation as ratio between travelable distance 
% to total original length
crushEfficiency = 0.6;  

crushResultStore = zeros(length(plateauForce), 3);

% crushable material's label
crushInd = 1;

% run simulations for each crushable material setting
for i = 1:length(plateauForce)

    % material properties for this crushable material design
    matProps = {plateauForce(i), mSeismometer};

    % select shock isolation model
    model = @crushModel;

    % get results from design point
    crushResultStore(i,1:2) = DP(data, model, matProps, debug, ...
        filterDebug, crushInd, masses, g, endTimes, energyThreshold, ...
        initialVsY, initialVsZ);

end

% apply crush efficiency
crushResultStore(:,1) = crushResultStore(:,1) / crushEfficiency;

% assign shock isolation label
crushResultStore(:,3) = ones(size(crushResultStore, 1), 1) * crushInd;

disp(['Finished crushable material models in ', ...
    num2str(toc(s)) , ' seconds'])

%% Assemble different paraemters to track between damper models

disp("Starting damper models")

s = tic;

% hydraulic damping constants [N/(m/s2)^2]
cVals = (10:100:1000).';

% space efficiency of shock isolation as ratio between travelable distance 
% to total original length
damperEfficiency = 0.5;

damperResultStore = zeros(length(cVals), 3);

% hydraulic damper's label
damperInd = 2;

dd = waitbar(0, ['Damper Model, elapsed time: ', num2str(toc(s)), ' seconds']);

for i = 1:length(cVals)

    % material properties for hydraulic damper (spring rate = 0 here)
    matProps = {0, cVals(i), mSeismometer};

    % shock isolation model
    model = @springDamperModel;

    % get results from this design point
    damperResultStore(i,1:2) = DP(data, model, matProps, debug, ...
        filterDebug, damperInd, masses, g, endTimes, energyThreshold, ...
        initialVsY, initialVsZ);

    % waitbar shows simulation progress
    waitbar(( i ) / length(cVals), dd, ...
            ['Damper Model, elapsed time: ', num2str(toc(s)), ' seconds'])
end

close(dd);

% apply space efficiency
damperResultStore(:,1) = damperResultStore(:,1) / damperEfficiency;

% add shock isolation label to data
damperResultStore(:,3) = ones(size(damperResultStore, 1), 1) * damperInd;

disp(['Finished damper models in ', num2str(toc(s)) , ' seconds'])

%% Spring damper model

disp("Starting spring damper models")

s = tic;

% hydraulic damping constants [N/(m/s2)^2]
cVals2 = (10:100:1000).';

% spring rates [N/m]
kVals2 = (1000:10000:100000).';

springDamperResultStore = zeros(length(cVals2), 3);

totalCombos = length(kVals2) * length(cVals2);

% spring hydraulic damper's label
springDamperInd = 3;

% at low damping values, the simulation time can be longer. using a waitbar
% to track
f = waitbar(0, ['Spring Damper Model, elapsed time: ', ...
    num2str(toc(s)), ' seconds']);

for j = 1:length(kVals2)
    for i = 1:length(cVals2)

        % spring damper material properties
        matProps = {kVals2(j), cVals2(i), mSeismometer};

        % select shock isolation model
        model = @springDamperModel;

        % store results from simulation
        springDamperResultStore((j-1) * length(cVals2) + i,1:2) = ...
            DP(data, model, matProps, debug, filterDebug, ...
            springDamperInd, masses, g, endTimes, energyThreshold, ...
            initialVsY, initialVsZ);

        % waitbar to track progress
        waitbar(((j-1) * length(cVals2) + i ) / totalCombos, f, ...
            ['Spring Damper Model (', num2str((j-1) * length(cVals2) + i), ...
            '/',num2str(totalCombos),') elapsed time: ', num2str(toc(s)), ' seconds'])
        
    end
end

close(f);

% apply space efficiency
springDamperResultStore(:,1) = springDamperResultStore(:,1) / damperEfficiency;

% apply shock isolation label to data
springDamperResultStore(:,3) = ones(size(springDamperResultStore, 1), 1) * springDamperInd;

disp(['Finished spring damper models in ', num2str(toc(s)) , ' seconds'])

%% get the pareto front input

% append results
paretoInput = [rigidResultStore; crushResultStore; damperResultStore; ...
    springDamperResultStore];

%% save pareto front data after going through all cases

save("paretoInput.mat", "paretoInput");

%% pareto front

load("paretoInput.mat");

% get the points along the Pareto front 
paretoFront = PF(paretoInput);

% number of design points for each shock isolation category
lRigid = 1;
lCrush = length(plateauForce);
lDamper = length(cVals);
lSpringDamper = length(cVals2) * length(kVals2);

% thesis plot display settings
markerSize = 300;
markerSizePareto = markerSize/3;
lineWidth = 1.5;
lineWidthLegend = 2;
markerSizeLegend = 40;
markerSizeLegend1 = 20;

% Figure 4-9 in Mike's thesis
f1 = figure;
f1.Position = [50,50, 700,500];
hold on
p1 = scatter(paretoInput(1:lRigid,1), paretoInput(1:lRigid,2) / g, ...
    markerSize, '*', 'LineWidth', lineWidth);
p2 = scatter(paretoInput(lRigid+1:lRigid+lCrush,1), ...
    paretoInput(lRigid+1:lRigid+lCrush,2) / g, markerSize, 'x', ...
    'LineWidth', lineWidth);
p3 = scatter(paretoInput(lRigid+lCrush+1:lRigid+lCrush+lDamper,1), ...
    paretoInput(lRigid+lCrush+1:lRigid+lCrush+lDamper,2) / g, ...
    markerSize, '+', 'LineWidth', lineWidth);
p4 = scatter(paretoInput(lRigid+lCrush+lDamper+1:lRigid+lCrush+lDamper+lSpringDamper,1), ...
    paretoInput(lRigid+lCrush+lDamper+1:lRigid+lCrush+lDamper+lSpringDamper,2) / g, markerSize, '.');
p5 = scatter(paretoFront(:,1), paretoFront(:,2) / g, markerSizePareto, 'o', 'LineWidth', lineWidth);
hold off
p1.MarkerEdgeColor = "#0072BD";
p2.MarkerEdgeColor = "#D95319";
p3.MarkerEdgeColor = "#77AC30";
p4.MarkerEdgeColor = "#7E2F8E";
p5.MarkerEdgeColor = "#000000";
xlabel("Shock Isolation Axial Length [m]")
ylabel("Max Axial Acceleration [g]")
% title("Shock Isolation Pareto Front")
[~, objh] = legend("Control","Crushable", "Damper", "Spring Damper", ...
    "Pareto Front", 'location','best', 'FontSize', 14);

% adjust legend element sizes
objhl = findobj(objh, 'type', 'patch'); 
set(objhl(4), 'MarkerSize', markerSizeLegend); 
set(objhl([1,2,3]), 'MarkerSize', markerSizeLegend1); 
set(objhl(5), 'MarkerSize', markerSizeLegend/3);
set(objhl, 'LineWidth', lineWidthLegend); 

ax = gca;
ax.FontSize = 14;
xlim([0 1.7])

%% show the pareto front with a limited view

% label the index of the desired design point
designPointInd = 9;

f2 = figure;
f2.Position = [750,50,700,500];
hold on
p1 = scatter(paretoInput(1:lRigid,1), paretoInput(1:lRigid,2) / g, ...
    markerSize, '*', 'LineWidth', lineWidth);
p2 = scatter(paretoInput(lRigid+1:lRigid+lCrush,1), ...
    paretoInput(lRigid+1:lRigid+lCrush,2) / g, markerSize, 'x', ...
    'LineWidth', lineWidth);
p3 = scatter(paretoInput(lRigid+lCrush+1:lRigid+lCrush+lDamper,1), ...
    paretoInput(lRigid+lCrush+1:lRigid+lCrush+lDamper,2) / g, ...
    markerSize, '+', 'LineWidth', lineWidth);
p4 = scatter(paretoInput(lRigid+lCrush+lDamper+1:lRigid+lCrush+lDamper+lSpringDamper,1), ...
    paretoInput(lRigid+lCrush+lDamper+1:lRigid+lCrush+lDamper+lSpringDamper,2) / g, markerSize, '.');
p5 = scatter(paretoFront(:,1), paretoFront(:,2) / g, ...
    markerSizePareto, 'o', 'LineWidth', lineWidth);
p6 = scatter(paretoFront(designPointInd,1), paretoFront(designPointInd,2) / g, markerSize, 'square', 'LineWidth', lineWidth);
hold off

p1.MarkerEdgeColor = "#0072BD";
p2.MarkerEdgeColor = "#D95319";
p3.MarkerEdgeColor = "#77AC30";
p4.MarkerEdgeColor = "#7E2F8E";
p5.MarkerEdgeColor = "#000000";
p6.MarkerEdgeColor = "#FF0000";

xlabel("Shock Isolation Axial Length [m]")
ylabel("Max Axial Acceleration [g]")
% title("Shock Isolation Pareto Front")

% legend properties
[a, objh] = legend("Control","Crushable", "Damper", "Spring Damper",...
    "Pareto Front", "Selected Design Point",'location','northwest', ...
    'FontSize', 14);
objhl = findobj(objh, 'type', 'patch'); 
set(objhl(4), 'MarkerSize', markerSizeLegend); 
set(objhl([1,2,3]), 'MarkerSize', markerSizeLegend1);
set(objhl(5), 'MarkerSize', markerSizeLegend/3); 
set(objhl(6), 'MarkerSize', markerSizeLegend*5/12); 
set(objhl, 'LineWidth', lineWidthLegend);

xlim([0 0.2])
ax1 = gca;
ax1.FontSize = 14;
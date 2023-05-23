%% Description

%{

    This design point (DP) function incorporates a given different material 
    design and calculates the responses from each impact case.

    Inputs:
    - data: cell array with all data from impact cases
    - model: material model's function
    - matProps: corresponding material properties for this design point
    - debug: 1 when plotting SRM data at each impact case
    - filterDebug: 1 when plotting how design point data is filtered to 
        one output result
    - modelInd: material model label
    - masses: 1x3 array of masses for each entity [kg]. This is in the
    format [body mass, flare mass, seismometer mass]    
    - g: gravity constant [m/s2]
    - endTimes: 1xN array with the impact durations [s] for each impact 
    case, where N is the number of impact cases
    - energyThreshold: threshold for when to stop simulation based on 
    seismometer energy [J]
    - initialVsY: 1xN array with the initial Y velocities [m/s] for each
    impact case, where N is the number of impact cases
    - initialVsZ: 1xN array with the initial Z velocities [m/s] for each
    impact case, where N is the number of impact cases

    Outputs:
    - designPointResult: 1x2 matrix, format = [maxStrokeDP, maxAccelDP]
        This is such that maxStrokeDP is the maximum stroke seen for
        this design point [m] and maxAccelDP is the maximum acceleration 
        of the seismometer for this design point [m/s2].

%}

%% Function

function designPointResult = DP(data, model, matProps, debug, ...
    filterDebug, modelInd, masses, g, endTimes, energyThreshold, ...
    initialVsY, initialVsZ)

    nCases = length(data);
    results = zeros(nCases, 4);

    mBody = masses(1);
    mFlare = masses(2);
    mSeismometer = masses(3);

    % run through all impact cases
    for i = 1:nCases
        
        % get the impact case data
        icData = data{i};
    
        % run the structural response model using the chosen material model
        [aPrime, vPrime, pPrime, breakInd, matEnergy] = SRM(model, ...
            matProps, icData, mBody, mFlare, mSeismometer, ...
            endTimes(i), initialVsY(i), initialVsZ(i));

        initMatEnergy = matEnergy;
        initYDisp = pPrime(end, 1) - pPrime(end, 3);
        initVsY = vPrime(end, 1) - vPrime(end, 3);
    
        % just define it...
        dt = 1e-4;

        % continue seismometer movement if needed after body stops moving
        [aPrime2, vPrime2, pPrime2, ~, ~] = SRMPost(model, matProps, ...
            mBody, mFlare, mSeismometer, initVsY, 0, initMatEnergy, ...
            initYDisp, energyThreshold, dt, modelInd);

        t = icData(1:breakInd,1);
        lenPost = size(aPrime2, 1);
        
        % check whether the seismometer is moving after the body stops
        if lenPost > 1
            tPost = dt * (0:lenPost-1) + t(end);
    
            % add the times from initial movement and post
            t = [t; tPost.'];
        
            % append entries from SRMPost if seismometer continues moving 
            % after body has stopped
            aPrime = [aPrime; aPrime2];
            vPrime = [vPrime; vPrime2];
            pPrime = [pPrime; pPrime2];

        end
       
        % get min/max accelerations of seismometer [m/s2]
        minAccel = min(aPrime(:, 1));
        maxAccel = max(aPrime(:, 1));
    
        % get min/max relative displacements[m]
        minRelDisp = min(pPrime(:, 1) - pPrime(:, 3));
        maxRelDisp = max(pPrime(:, 1) - pPrime(:, 3));
    
        % generate result vector for each impact case: extrema for seismometer
        % axial acceleration and extrema for relative displacements between
        % body and seismometer (negative indicates compression)
        results(i, :) = [minAccel, maxAccel, minRelDisp, maxRelDisp];
        
        % debugging
        if debug
            figure
            subplot(2,2,1)
            hold on
            plot(t, aPrime(:,3) / g)
            plot(t, aPrime(:,1) / g)
            hold off
            xlabel('time [s]')
            ylabel('axial acceleration [g]')
            legend('penetrator body', 'seismometer', 'Location','best')
            
            subplot(2,2,2)
            hold on
            plot(t, (pPrime(:, 3) - pPrime(:, 1)) * 100)
            hold off
            xlabel('time [s]')
            ylabel('axial movement [cm]')
            
            subplot(2,2,3)
            hold on
            plot(t, (vPrime(:, 3) - vPrime(:, 1)) )
            hold off
            xlabel('time [s]')
            ylabel('axial velocity [m/s]')
            
            subplot(2,2,4)
            hold on
            plot(t, (vPrime(:, 4) - vPrime(:, 2)) )
            hold off
            xlabel('time [s]')
            ylabel('lateral velocity [m/s]')
        
            sgtitle("Impact case " + i)
        end
    end
       
    % if multiple points with same max acceleration, take the one with max stroke
    [maxAccelDP, maxAccelInd] = max(results(:,2));
    [maxStrokeDP, maxStrokeInd] = max(results(:,4) - results(:,3));

    designPointResult = [maxStrokeDP, maxAccelDP];
    
    % debuggings
    if filterDebug

        % scale dx to be appropriate
        xMax = max(results(:,4) - results(:,3));
        xMin = min(results(:,4) - results(:,3));
    
        % set label position
        dx = (xMax - xMin)/ 50;
        dy = 0;
        
        % filtering debug plot
        figure
        hold on
        scatter(results(:,4) - results(:,3), results(:,2) / g)
        text(results(:,4) - results(:,3) + dx, results(:,2)/g + dy, ...
            cellstr(num2str((1:nCases).')), 'FontSize', 12)
        scatter(maxStrokeDP, results(maxStrokeInd,2)/g)
        scatter(results(maxAccelInd,4) - results(maxAccelInd,3), ...
            maxAccelDP / g)
        hold off
        xlabel("Stroke Movement [m]")
        ylabel("Max Acceleration [g]")
        title("Design point's max acceleration and stroke distance for each impact case")
        legend("all impact cases", "max stroke case", "max accel case", ... 
            'location','best')
    end
end
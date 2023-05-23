%% Description

%{

    Structural response model function for after the body stops. It is
    assumed for simplicity that all accelerations and motion are along the
    Y' axis and the body remains fixed. Thus, the simulation is run
    assuming no body accelerations and assuming 0 degree tipping angle.

    Firn coordinate system = FCS, body coordinate system = BCS.

    Inputs:
    - model: shock isolation model function
    - matProps: properties needed for shock isolation model
    - mBody: body mass [kg]
    - mFlare: flare mass [kg]
    - mSeismometer: seismometer mass [kg]
    - endTime: impact duration [s]
    - initialVsY: initial Y velocity in this post impact region [m/s]
    - initialVsZ: initial Z velocity in this post impact region [m/s]
    - initMatEnergy: energy stored in shock isolation [J] at start of this 
    post impact region
    - initYDisp: relative displacement between seismometer and body [m] at
    start of this post impact region
    - energyThreshold: seismometer energy treshold for deciding when to end
    simulation [J]
    - dt: timestep size [s]
    - modelInd: shock isolation model label

    Outputs:
    - aPrime: nx4 matrix storing accelerations of body and seismometer in
    BCS [m/s2]. Format is [seismometer Y' accel, seismometer Z' accel, 
    body Y' accel, body Z' accel]
    - vPrime: nx4 matrix storing velocities of body and seismometer in
    BCS [m/s]. Format is [seismometer Y' velocity, seismometer Z' velocity, 
    body Y' velocity, body Z' velocity]
    - pPrime: nx4 matrix storing positions of body and seismometer in
    BCS [m]. Format is [seismometer Y' position, seismometer Z' position, 
    body Y' position, body Z' position]
    - breakInd: index at which the simulation loop stopped
    - matEnergy: energy stored in shock isolation [J]

%}

%% Function 

function [aPrime, vPrime, pPrime, breakInd, matEnergy] = SRMPost(model, matProps, mBody, ...
    mFlare, mSeismometer, initialVsY, initialVsZ, ...
    initMatEnergy, initYDisp, energyThreshold, dt, modelInd)
    
    % zero body acceleration (assumption)
    accZ = 0;
    accY = 0;
    
    %% Run through states
    
    % for storing the state values during the integration
    stateStore = zeros(1, 15);
    
    mTot = mBody + mFlare;
    
    % assume now the penetrator has zero starting velocity, and the
    % seismometer is traveling at its ending axial velocity in the Y
    % direction
    % state in the format: [ys, ysdot, zs, zsdot, yp, ypdot, zp, zpdot];
    statei = [initYDisp, initialVsY, 0, 0, 0, ...
        0, 0, 0];

    breakInd = 1;

    % energy stored in shock isolation
    matEnergy = initMatEnergy;
    
    i=0;

    % booleans that trigger when max and min positions are found
    maxFind = 0;
    minFind = 0;

    % continue simulating until end conditions are met
    while true
        i = i+1;

        % kinetic energy of seismometer and energy stored in shock
        % isolation
        seismometerEnergy = 1/2*mSeismometer*statei(2)^2 + matEnergy;

        % stop simulation if seismometer energy is below a threshold
        if seismometerEnergy <= energyThreshold
            % store the index before for the purpose of plotting only the
            % region of interest
            breakInd = i - 1;
    
            break;
        end
        
        % assume for now that the Y' axis is aligned with the Y axis
        theta = 0;
    
        [ayprime, azprime] = YZtoYZPrime(theta, accY, accZ);

        % get the relative axial velocity between the body and seismometer in
        % the prime frame
        relVelYPrime = (statei(2) * cosd(theta) + statei(4) * sind(theta)) - ...
            (statei(6) * cosd(theta) + statei(8) * sind(theta));

        % get the relative axial position between the body and seismometer in
        % the prime frame
        relPosYPrime = (statei(1) * cosd(theta) + statei(3) * sind(theta)) - ...
            (statei(5) * cosd(theta) + statei(7) * sind(theta));

        % shock isolation model
        [ayprimes, matEnergy] = model(ayprime, relVelYPrime, relPosYPrime, ...
            matProps);

        [ays, azs, ayb, azb] = RBMM(mTot, mSeismometer, ayprime, ...
        ayprimes, azprime, theta);

        %%%% integrator is below

        % use forward euler to move the next state (only along Y' direction
        % for the seismometer)
        statei1 = statei + dt * [statei(2), ays, 0, 0, ...
            0, 0, 0, 0];
   
        % the above is the predicted next state. now, adjust using the
        % adjusting stage of Heun's method
        vi1Predicted = statei1(2:2:end);
        xi = statei(1:2:end-1);

        % the adjusted position
        xi1 = xi + dt / 2 * (statei(2:2:end) + vi1Predicted);
    
        % now, calculate the predicted acceleration given the predicted
        % state (from statei1)
            
        thetai1 = 0;

        [ayprimei1, azprimei1] = YZtoYZPrime(thetai1, accY, accZ);

        % get the relative axial velocity between the body and seismometer in
        % the prime frame
        relVelYPrimei11 = (statei1(2) * cosd(thetai1) + statei1(4) * sind(thetai1)) - ...
            (statei1(6) * cosd(thetai1) + statei1(8) * sind(thetai1));

        % get the relative axial position between the body and seismometer in
        % the prime frame
        relPosYPrimei11 = (statei1(1) * cosd(thetai1) + statei1(3) * sind(thetai1)) - ...
            (statei1(5) * cosd(thetai1) + statei1(7) * sind(thetai1));

        [ayprimesi1, ~] = model(ayprimei1, relVelYPrimei11, relPosYPrimei11, ...
            matProps);

        [aysi1, ~, ~, ~] = RBMM(mTot, mSeismometer, ayprimei1, ...
        ayprimesi1, azprimei1, thetai1);

        % the accelerations of seismometer and body in XY coordinates
        accelsi = [ays, 0, 0, 0];
        accelsi1P = [aysi1, 0, 0, 0];

        vi = statei(2:2:end);

        % the adjusted velocity
        vi1 = vi + dt / 2 * (accelsi + accelsi1P);

        statei1(1:2:end-1) = xi1;
        statei1(2:2:end) = vi1;

        %%%% integrator is above

        thetai1 = 0;

        relVelYPrimei1 = (statei1(2) * cosd(thetai1) + statei1(4) * sind(thetai1)) - ...
                (statei1(6) * cosd(thetai1) + statei1(8) * sind(thetai1));

        
        % store first state if i == 1
        if i == 1
            stateStore = [statei, ays, azs, ayb, azb, ...
                 relVelYPrime, relVelYPrimei1, theta];
        else
            stateStore = [stateStore; statei, ays, azs, ayb, azb, ...
                 relVelYPrime, relVelYPrimei1, theta];
        end

        % only track these if the model is a spring damper
        if modelInd == 3 || modelInd == 4
            
            % find the maximum compression
            if statei1(2) > 0 && statei(2) <= 0

                maxFind = 1;

            % find the minimum compression (max extension)
            elseif statei1(2) < 0 && statei(2) >= 0
               
                minFind = 1;

            end
        
            % break from the integration loop if so
            if maxFind && minFind

                breakInd = i - 1;
   
                break;
            end 
            
        end

        % reset the state
        statei = statei1;
    
    end

    %% convert to prime coordinates
    
    aSeisArr = [stateStore(:, 9).'; stateStore(:, 10).'];
    aBodyArr = [stateStore(:, 11).'; stateStore(:, 12).'];
    
    vSeisArr = [stateStore(:, 2).'; stateStore(:, 4).'];
    vBodyArr = [stateStore(:, 6).'; stateStore(:, 8).'];
    
    pSeisArr = [stateStore(:, 1).'; stateStore(:, 3).'];
    pBodyArr = [stateStore(:, 5).'; stateStore(:, 7).'];
    
    aPrime = zeros(size(aBodyArr,2), 4);
    vPrime = zeros(size(vBodyArr,2), 4);
    pPrime = zeros(size(pBodyArr,2), 4);
    
    % for i = 1:length(aPrime)
    for i = 1:size(aPrime, 1)
    
        theta = stateStore(i, end);
    
        Rrot = [cosd(theta), sind(theta); -sind(theta), cosd(theta)];
    
        % convert acceleration
        aSeisPrime = Rrot * aSeisArr(:,i);
        aBodyPrime = Rrot * aBodyArr(:,i);
    
        aPrime(i,:) = [aSeisPrime(1), aSeisPrime(2), aBodyPrime(1), aBodyPrime(2)];
    
        % convert velocity
        vSeisPrime = Rrot * vSeisArr(:,i);
        vBodyPrime = Rrot * vBodyArr(:,i);
    
        vPrime(i, :) = [vSeisPrime(1), vSeisPrime(2), vBodyPrime(1), vBodyPrime(2)];
    
        % convert position
        pSeisPrime = Rrot * pSeisArr(:,i);
        pBodyPrime = Rrot * pBodyArr(:,i);
    
        pPrime(i,:) = [pSeisPrime(1), pSeisPrime(2), pBodyPrime(1), pBodyPrime(2)];
    
    end

end
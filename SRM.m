%% Description

%{

    Structural response model function. Firn coordinate system = FCS, body
    coordinate system = BCS.

    Inputs:
    - model: shock isolation model function
    - matProps: properties needed for shock isolation model
    - data: impact case data
    - mBody: body mass [kg]
    - mFlare: flare mass [kg]
    - mSeismometer: seismometer mass [kg]
    - endTime: impact duration [s]
    - initialVsY: initial Y velocity [m/s]
    - initialVsZ: initial Z velocity [m/s]

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

function [aPrime, vPrime, pPrime, breakInd, matEnergy] = SRM(model, ...
    matProps, data, mBody, mFlare, mSeismometer, endTime, initialVsY, ...
    initialVsZ)
    
    accZ = data(:,6);
    accY = data(:,2);
    thetaInt = data(:,8);
    
    %% Run through states
    
    % time [s]
    t = data(:,1);
    
    % for storing the state values during the integration
    stateStore = zeros(length(data(:,1)), 15);
    
    % total mass [kg]
    mTot = mBody + mFlare;
    
    % get the change in time for each timestep (assume centered around the
    % point in time
    dt = diff(t);
    
    % append an extra dt at the end so vectors are consistent
    dt = [dt; dt(end)];
    
    % state in the format: [ys, ysdot, zs, zsdot, yp, ypdot, zp, zpdot];
    statei = [0, initialVsY, 0, initialVsZ, 0, ...
        initialVsY, 0, initialVsZ];
    
    breakInd = 1;
    
    % loop through all timesteps
    for i = 1:length(t)

        % break from current loop when penetrator has reached peak depth
        if t(i+1) >= endTime
            % store the index before for the purpose of plotting only the
            % region of interest
            breakInd = i - 1;
    
            break;
        end
        
        % current timestep's penetrator tipping angle [deg]
        theta = thetaInt(i);
    
        % convert body penetration accelerations from FCS to BCS
        [ayprime, azprime] = YZtoYZPrime(theta, accY(i), accZ(i));

        % get the relative axial velocity between the body and seismometer in
        % the prime frame
        relVelYPrime = (statei(2) * cosd(theta) + statei(4) * sind(theta)) - ...
            (statei(6) * cosd(theta) + statei(8) * sind(theta));

        % get the relative axial position between the body and seismometer in
        % the prime frame
        relPosYPrime = (statei(1) * cosd(theta) + statei(3) * sind(theta)) - ...
            (statei(5) * cosd(theta) + statei(7) * sind(theta));

        [ayprimes, ~] = model(ayprime, relVelYPrime, relPosYPrime, ...
            matProps);

        [ays, azs, ayb, azb] = RBMM(mTot, mSeismometer, ayprime, ...
        ayprimes, azprime, theta);
        
        %%%% integrator is below

        % use forward euler to move to the next state
        statei1 = statei + dt(i) * [statei(2), ays, statei(4), azs, ...
            statei(6), ayb, statei(8), azb];

        % the above is the predicted next state. now, adjust using the
        % adjusting stage of Heun's method
        vi1Predicted = statei1(2:2:end);
        xi = statei(1:2:end-1);

        % the adjusted position
        xi1 = xi + dt(i) / 2 * (statei(2:2:end) + vi1Predicted);
    
        % now, calculate the predicted acceleration given the predicted
        % state (from statei1)
            
        thetai1 = thetaInt(i+1);

        [ayprimei1, azprimei1] = YZtoYZPrime(thetai1, accY(i+1), accZ(i+1));

        % get the relative axial velocity between the body and seismometer in
        % the prime frame
        relVelYPrimei11 = (statei1(2) * cosd(thetai1) + statei1(4) * sind(thetai1)) - ...
            (statei1(6) * cosd(thetai1) + statei1(8) * sind(thetai1));

        % get the relative axial position between the body and seismometer in
        % the prime frame
        relPosYPrimei11 = (statei1(1) * cosd(thetai1) + statei1(3) * sind(thetai1)) - ...
            (statei1(5) * cosd(thetai1) + statei1(7) * sind(thetai1));

        
        [ayprimesi1, matEnergy] = model(ayprimei1, relVelYPrimei11, relPosYPrimei11, ...
            matProps);

        [aysi1, azsi1, aybi1, azbi1] = RBMM(mTot, mSeismometer, ayprimei1, ...
        ayprimesi1, azprimei1, thetai1);

        % the accelerations of seismometer and body in XY coordinates
        accelsi = [ays, azs, ayb, azb];
        accelsi1P = [aysi1, azsi1, aybi1, azbi1];

        vi = statei(2:2:end);

        % the adjusted velocity
        vi1 = vi + dt(i) / 2 * (accelsi + accelsi1P);

        statei1(1:2:end-1) = xi1;
        statei1(2:2:end) = vi1;

        %%%% integrator is above    

        thetai1 = thetaInt(i+1);

        % relative velocity between body and seismometer at next timestep
        relVelYPrimei1 = (statei1(2) * cosd(thetai1) + statei1(4) * sind(thetai1)) - ...
                (statei1(6) * cosd(thetai1) + statei1(8) * sind(thetai1));

        stateStore(i, :) = [statei, ays, azs, ayb, azb, ...
        relVelYPrime, relVelYPrimei1, theta];

        % reset the state
        statei = statei1;
    
    end

    % eliminate the sections of stateStore that aren't used
    stateStore(breakInd+1:end, :) = [];
    
    %% Convert variables to prime coordinates
    
    aSeisArr = [stateStore(:, 9).'; stateStore(:, 10).'];
    aBodyArr = [stateStore(:, 11).'; stateStore(:, 12).'];
    
    vSeisArr = [stateStore(:, 2).'; stateStore(:, 4).'];
    vBodyArr = [stateStore(:, 6).'; stateStore(:, 8).'];
    
    pSeisArr = [stateStore(:, 1).'; stateStore(:, 3).'];
    pBodyArr = [stateStore(:, 5).'; stateStore(:, 7).'];
    
    aPrime = zeros(size(aBodyArr,2), 4);
    vPrime = zeros(size(vBodyArr,2), 4);
    pPrime = zeros(size(pBodyArr,2), 4);
    
    for i = 1:length(aPrime)
    
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
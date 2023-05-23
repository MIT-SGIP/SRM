%% Description

%{

    The rigid body motion model (RBMM) used to describe the accelerations
    of the seismometer and the body in the (Y', Z') coordinate system.

    Inputs:
    - mTot: total mass of body, flare, and seismometer [kg]
    - mSeismometer: mass of seismometer [kg]
    - ayprime: Y' acceleration from Ansys impact simulation [m/s2]
    - ayprimes: seismometer Y' acceleration [m/s2]
    - azprime: body and seismometer Z' acceleration [m/s2]
    - theta: penetrator tipping angle [deg]

    Outputs:
    - ays: seismometer Y acceleration [m/s2]
    - azs: seismometer Z acceleration [m/s2]
    - ayb: body Y acceleration [m/s2]
    - azb: body Z acceleration [m/s2]

%}

%% Function

function [ays, azs, ayb, azb] = RBMM(mTot, mSeismometer, ayprime, ...
        ayprimes, azprime, theta)    

        % rotation matrixs
        Rrot = [cosd(theta), sind(theta); -sind(theta), cosd(theta)];
    
        % net Y' acceleration of body [m/s2]
        netAyp = (mTot * ayprime -  mSeismometer * ayprimes) / ...
            (mTot - mSeismometer); 
        
        % group BCS accelerations by seismometer and body
        aPrimeSeis = [ayprimes; azprime];
        aPrimeBody = [netAyp; azprime];

        % transform coordinates back to the unprimed frame
        aSeis = Rrot\aPrimeSeis;
        aBody = Rrot\aPrimeBody;

        % split FCS acceleration outputs into separate variables
        ays = aSeis(1);
        azs = aSeis(2);
        ayb = aBody(1);
        azb = aBody(2);

end
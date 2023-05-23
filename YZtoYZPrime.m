%% Description

%{

    Rotate the accelerations from the firn coordinate system to the body
    coordinate system (prime coordinates).
    
    Inputs:
    - theta: tipping angle of body [deg]
    - accYi: Y acceleration [m/s2]
    - accZi: Z acceleration [m/s2]  

    Outputs:
    - ayprime: Y' acceleration [m/s2]
    - azprime: Z' acceleration [m/s2]  

%}

%% Function

function [ayprime, azprime] = YZtoYZPrime(theta, accYi, accZi)

        % rotation matrix
        Rrot = [cosd(theta), sind(theta); -sind(theta), cosd(theta)];
    
        % transform the accelerations into the penetrator body frame
        a = [accYi; accZi];

        aPrime = Rrot * a;

        % split data into components
        ayprime = aPrime(1);
        azprime = aPrime(2);

end
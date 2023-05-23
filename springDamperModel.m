%% Description

%{

    Spring damper model for shock isolation study. In this model, the spring
    and damper act in parallel. The spring has a linear spring rate, and
    the damper's force is related to the square of the piston's velocity
    due to the orifice effect.

    Input:
    - relVelYPrime: Y' relative velocity between the seismometer and body
    [m/s]. Negative when damper is compressing.
    - relPosYPrime: Y' relative position between the seismometer and body
    [m/s]. Negative when damper is compressing.
    - inputs: cell with three relevant properties for the material model
    with format {spring rate [N/m], damping coefficient [N/(m/s)^2], 
    seismometer mass [kg]}

    Output:
    - ayprimes: seismometer Y' acceleration [m/s2]
    - storedEnergy: energy stored in shock isolation but not
    dissipated [J]. This is the energy stored in the spring.
    
%}

%% Function

function [ayprimes, storedEnergy] = springDamperModel(~, relVelYPrime, ...
    relPosYPrime, inputs)
       
        % load material properties into local variables
        k = inputs{1};
        c = inputs{2};
        mSeismometer = inputs{3};

        % damping force [N]
        Fdamp = -c * sign(relVelYPrime) * relVelYPrime^2;

        % spring force [N]
        Fspring = -k * relPosYPrime;

        % spring and damper act in parallel
        Ftotal = Fspring + Fdamp;

        ayprimes = Ftotal / mSeismometer;

        % energy stored in spring but not yet dissipated [J]
        storedEnergy = 1/2 * k * (relPosYPrime).^2;
        
end
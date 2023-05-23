%% Description

%{

    Crushable material model for shock isolation study. In this model, the 
    crushable material is modelled as a rigid perfectly plastic material.

    Input:
    - ayprime: Y' acceleration of the body due to impact [m/s2]
    - relVelYPrime: Y' relative velocity between the seismometer and body
    [m/s]. Negative when crushable is compressing.
    - inputs: cell with twp relevant properties for the material model
    with format {crushable plateau force [N], seismometer mass [kg]}

    Output:
    - ayprimes: seismometer Y' acceleration [m/s2]
    - storedEnergy: energy stored in shock isolation but not
    dissipated [J]. Zero for crushable model used (no elastic deformation)
    
%}

%% Function

function [ayprimes, storedEnergy] = crushModel(ayprime, relVelYPrime, ...
    ~, inputs)
       
        % load material properties into local variables
        plateauForce = inputs{1};
        mSeismometer = inputs{2};
        
        % crush if the two bodies are approaching eachother or the
        % seismometer's inertial weight reaches the plateau force
        if relVelYPrime < 0 || ayprime > plateauForce / mSeismometer

            ayprimes = plateauForce / mSeismometer;
            
        % else, just use the body's Y' acceleration
        else

            ayprimes = ayprime;

        end

        % crushable model does not store potential energy
        storedEnergy = 0;

end
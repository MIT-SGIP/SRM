%% Description

%{

    Control for shock isolation study. In this model, the seismometer
    experiences the same Y' acceleration as the body.

    Input:
    - ayprime: Y' acceleration of the body due to impact [m/s2]

    Output:
    - ayprimes: seismometer Y' acceleration [m/s2]
    - storedEnergy: energy stored in shock isolation but not
    dissipated [J]. Always zero for control
    
%}

%% Function

function [ayprimes, storedEnergy] = controlModel(ayprime, ~, ~, ~)
        
        ayprimes = ayprime;

        storedEnergy = 0;
       
end
classdef Spacecraft
    properties
        state (6,:) double   % [position (km), velocity (km/s)]
        epoch (1,:) double   % Julian date or POSIX
        mass                 % wet mass (kg)
        Ad                   % wetted area, drag (m^2)
        As                   % absorbing area, SRP (m^2)
        Cd                   % drag coefficient
        Cr                   % reflectivity coefficient
        CB string            % central body
    end


    methods
        % Constructor
        function SC = Spacecraft(state,epoch,CB)
            SC.state = state;
            SC.epoch = epoch;
            SC.CB = Body(CB,epoch);
            
        end

    end






end
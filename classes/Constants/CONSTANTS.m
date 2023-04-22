classdef CONSTANTS
    properties (Constant)
        AU = 149597870.700 %[km] astronomical unit in km
        c = 299792.458; %[km/s] speed of light
        JD = 86400; %[s] Julian day in seconds
        JY = 365.25; %[d] Julian year in days
        JC = 36525; %[d] Julian century in days
        JD2000 = 2451545.0; % [JD] standard epoch of J2000.0
        GMconv = CONSTANTS.AU^-3*CONSTANTS.JD^2; % [AU^3/day^2]^-1 conversion of gravitational parameter between metric and astronomical units
        G = 6.67430e-20*CONSTANTS.GMconv; % [kg*AU^3/day^2] gravitational constant
        mu_sun = 1.32712440041279419e11*CONSTANTS.GMconv; % [AU^3/day^2] heliocentric gravitational constant
    end



end
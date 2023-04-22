classdef Body < handle
    properties
        name
        epoch %[JD] J2000 Julian date
    end

    properties (Dependent)
        mu
        kep struct
        equRadius
        meanRadius
        parent
        period
        meanMotion
        specificEnergy
        SOI
    end

    properties (Access = private)
        % coefficients for polynomial fit to orbital elements (J2000)
        % rows: [ meanLon,a,e,i,RAAN,lonPer]
        % cols: [ T0,T1,T2,T3 ]
        % source: "Meeus, Jean. 1991. Astronomical Algorithms. P.212-5."
        polyT (6,4) double 
    end

    
    methods
    % Constructor method
    function obj = Body(name,epoch)
            obj.name = name;
            obj.epoch = epoch;
            GMconv = CONSTANTS.GMconv;


            switch name
                case 'Mercury'
                    obj.mu = 22031.868551*GMconv; %[AU^3/d^2]
                    obj.parent = 'Sun';
                    obj.polyT = [
                        252.250906 1.494726746358e5 -5.36e-6 2e-9;  % L (deg)
                        0.387098310 0 0 0;                          % a (AU)
                        0.20563175 2.0407e-5 -2.83e-8 1.8e-10;      % e
                        7.004986 -5.9516e-3 8e-7 4.3e-8;            % i (deg from ecliptic)
                        48.330893 -0.1254227 -8.833e-5 -2e-7        % RAAN (deg)
                        77.456119 0.1588643 -1.342e-5 7e-9;         % lonPer (deg)
                    ];


                case 'Venus'
                    obj.mu = 324858.592000*GMconv; %[AU^3/d^2]
                    obj.parent = 'Sun';
                    obj.polyT = [
                        181.979801 5.85178156760e4 1.65e-6 -2e-9;   % L (deg)
                        0.723329820 0 0 0;                          % a (AU)
                        6.77192e-3 -4.7765e-4 9.81e-8 4.6e-10;      % e
                        3.394662 -8.568e-4 -3.244e-5 9e-9;          % i (deg from ecliptic)
                        76.679920 -0.2780134 -1.4257e-4 -1.64e-7    % RAAN (deg)
                        131.563703 4.8746e-3 -1.38467e-3 -5.695e-6; % lonPer (deg)
                    ];

                case 'Earth'
                    obj.mu = 398600.435507*GMconv; %[AU^3/d^2]
                    obj.parent = 'Sun';
                    obj.polyT = [
                        100.466457 3.59993728565e4 -5.68e-6 -1e-9;  % L (deg)
                        1.000001018 0 0 0;                          % a (AU)
                        1.670863e-2 -4.2037e-5 -1.267e-7 1.4e-9;    % e
                        0 1.30548e-2 -9.31e-6 3.4e-9;               % i (deg from ecliptic)
                        174.873176 -0.2410908 4.262e-5 1e9;         % RAAN (deg)
                        102.937348  0.3225654 1.4799e-4 -3.9e-8     % lonPer (deg)
                    ];

                case 'Moon'
                    obj.mu = 4902.800118*GMconv; %[AU^3/d^2]
                    obj.parent = 'Earth';

                case 'Mars'
                    obj.mu = 42828.375816*GMconv; %[AU^3/d^2]
                    obj.parent = 'Sun';

                case 'Jupiter'
                    obj.mu = 126712764.100000*GMconv; %[AU^3/d^2]
                    obj.parent = 'Sun';

                case 'Saturn'
                    obj.mu =  37940584.841800*GMconv; %[AU^3/d^2]
                    obj.parent = 'Sun';

                case 'Uranus'
                    obj.mu = 5794556.400000*GMconv; %[AU^3/d^2]
                    obj.parent = 'Sun';

                case 'Neptune'
                    obj.mu = 6836527.100580*GMconv; %[AU^3/d^2]
                    obj.parent = 'Sun';

                case 'Pluto'
                    obj.mu = 975.500000*GMconv; %[AU^3/d^2]
                    obj.parent = 'Sun';

            end
            
        end


    % Calculate dependent properties
        function Tp = get.period(obj)
            obj.kep = get.kep(obj);
            Tp = 2*pi*sqrt(obj.kep.a^3/obj.mu);

        end

        function n = get.meanMotion(obj)
            obj.kep = get.kep(obj);
            n = sqrt(obj.mu/obj.kep.a^3);
        end
        
        function E = get.specificEnergy(obj)
            obj.kep = get.kep(obj);
            E = -0.5*obj.mu/obj.kep.a;
        end

        function rSOI = get.SOI(obj)
            obj.kep = get.kep(obj);
            rSOI = obj.kep.a*(obj.mu/obj.parent.mu);
        end

        function T = calcTime(obj)
            % time measured in julian centuries from epoch of J2000 
            T = (obj.epoch - CONSTANTS.JD2000)/CONSTANTS.JC;
        end
        
        function kep = get.kep(obj)

        end
    
    end




end
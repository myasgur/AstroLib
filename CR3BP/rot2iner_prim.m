function [Xm1i,Xm2i] = rot2iner_prim(t,M,mu,del,spatial)

% Xin=rot2iner_prim(X,T,M,mu,del,spatial);
%
% Transformation from M body-centered rotating (nondim.) coordinates to
% inertial (nondim.) coordinated
%
% M = 2 : smaller mass(M2) centered inertial coordinates
% M = 1 : LARGER  mass(M1) centered inertial coordinates
% M = 0 : center-of-mass   centered inertial coordinates
%
% t   = nondimensional CR3BP time
% del = the rotation offset at time t=0 (in radians)
% spatial = (boolean) determine the size of the state (4DOF vs 6DOF)
%
% Xrot = inv(B)*(Xin-A)
% see p.6 of Cassall(1996)
%
% NOTE: nondim. means nondimensional CR3BP units where
%   sum of primaries' mass = 1;
%   constant distance between primaries = 1;
%   period of primaries' orbit = 2*pi
%
%-------------------------------------------------------------------------
% CR3BP (Circular Restricted Three-Body [Gravitational] Problem)
% with the LARGER MASS, M1, to the left of the origin at (-mu,0)
% and the smaller mass, M2, or the planet (i.e., Earth) at (1-mu,0)
%
%   
%       (rotating coords)
%
%                L4
% -L3----m1--+-----L1--m2--L2-
%                L5
%

% rotating positions of the two primaries
Xm1rot = repmat([ -mu, 0],size(t));
Xm2rot = repmat([1-mu, 0],size(t));

% Defaults
if nargin<=5
    spatial=0;
    if nargin<=4
        del=0;
        if nargin<=3
            mu=0;
            if nargin<=2
                M=0;
            end
        end
    end
end

if     M==1, d= -mu;        % distance to LARGER  primary in CR3BP (M1)
elseif M==2, d=1-mu;        % distance to smaller primary in CR3BP (M2)
elseif M==0, d=0;           % center-of-mass is the origin
end

c=cos(t + del); s = sin(t + del);

Xm1i = [];
Xm1i(:,1)= c.*(Xm1rot(:,1)-d) - s.*Xm1rot(:,2);
Xm1i(:,2)= s.*(Xm1rot(:,1)-d) + c.*Xm1rot(:,2);

Xm2i = [];
Xm2i(:,1)= c.*(Xm2rot(:,1)-d) - s.*Xm2rot(:,2);
Xm2i(:,2)= s.*(Xm2rot(:,1)-d) + c.*Xm2rot(:,2);

if spatial
    % rotating about z axis
    Xm1i(:,3) = zeros(size(t));
    Xm2i(:,3) = zeros(size(t));
end

end
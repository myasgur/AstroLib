function Xin = rot2iner(X,t,M,mu,del)
% Xin=rot2iner(X,T,M,mu,del);
%
% Transformation from M body-centered rotating (nondim.) coordinates to
% inertial (nondim.) coordinated
%
% M = 2 : smaller mass(M2) centered inertial coordinates
% M = 1 : LARGER  mass(M1) centered inertial coordinates
% M = 0 : center-of-mass   centered inertial coordinates
%
% X   = state 4- or 6-vector matrix (in nondim. rotating coordinates)
% t   = nondimensional CR3BP time
% del = the rotation offset at time t=0 (in radians)
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

% Defaults
if nargin<=4
    del=0;
    if nargin<=3
        mu=0;
        if nargin<=2
            M=0;
        end
    end
end

if     M==1, d= -mu;        % distance to LARGER  primary in CR3BP (M1)
elseif M==2, d=1-mu;        % distance to smaller primary in CR3BP (M2)
elseif M==0, d=0;           % center-of-mass is the origin
end

c=cos(t + del); s = sin(t + del);
x = 1; y = 2;               % indices 

Xin = [];
Xin(:,x)= c.*(X(:,x)-d) - s.*X(:,y);
Xin(:,y)= s.*(X(:,x)-d) + c.*X(:,y);

if size(X,2)==4
    vx = 3; vy = 4;
elseif size(X,2)==6
    vx = 4; vy = 5;

    % rotating about z axis
    Xin(:,3) = X(:,3);
    Xin(:,6) = X(:,6);
end

Xin(:,vx) =-s.*(X(:,x)+X(:,vy)-d) + c.*(X(:,vx) - X(:,y));
Xin(:,vy) = c.*(X(:,x)+X(:,vy)-d) + s.*(X(:,vx) - X(:,y));

end
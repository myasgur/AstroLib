function E = energy(x,mu)
%
%        E = energy(x,mu)
%
% with the LARGER MASS, m1 to the left of the origin at (-mu,0)
% and m2, or the planet (ie. Earth), is at (1 - mu, 0)
%
%                L4
% -L3----m1--+-----L1--m2--L2-
%                L5
%
mu1 = 1-mu; % mass of larger  primary (nearest origin on left)
mu2 =   mu; % mass of smaller primary (furthest from origin on right)

dim = size(x,2);


switch dim
    case 4
        r= ((x(:,1)+mu2).^2 + x(:,2).^2).^0.5;  % r: distance to m1, LARGER MASS
        R= ((x(:,1)-mu1).^2 + x(:,2).^2).^0.5;  % R: distance to m2, smaller mass

        % Kinetic Energy
        Vsq = (x(:,3).^2 + x(:,4).^2);
    case 6
        r= ((x(:,1)+mu2).^2 + x(:,2).^2 + x(:,3).^2).^0.5; 
        R= ((x(:,1)-mu1).^2 + x(:,2).^2 + x(:,3).^2).^0.5; 

        % Kinetic Energy
        Vsq = (x(:,4).^2 + x(:,5).^2 + x(:,6).^2);
end

% Psuedo-potential
U = -0.5*(x(:,1).^2 + x(:,2).^2) - (mu1./r + mu2./R) - 0.5*mu1*mu2;

% Jacobi constant (total system energy)
E = 0.5*Vsq + U;

end
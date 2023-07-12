function [L,gamma,xL] = findGamma(Lpt,mu)

if nargin < 2
    mu = 1/82.3;
    if nargin < 1
        Lpt = 1;
    end
end

tol = eps(1);

%preallocate arrays for speed
N = length(Lpt);
L = zeros(1,N);
gamma = zeros(1,N);
xL = zeros(1,N);

% solve for x-position of collinear Lagrange points using iterative scheme
for i = 1:N
    if Lpt(i) == 1 || Lpt(i) == 2
        gam0 = (mu*(1-mu)/3)^(1/3);
        gam = gam0 + 1;
        p = (-1)^Lpt(i); % -1 (L1) or +1 (L2)

        while abs(gam - gam0) > tol
            gam0 = gam;
            gam = (mu*(gam0+p)^2 / (3-2*mu+p*gam0*(3-mu+p*gam0)))^(1/3);
        end
        gamma(i) = gam;
        xL(i) = 1 - mu + p*gam;

    elseif Lpt(i) == 3

        gam0 = (mu*(1-mu)/3)^(1/3);
        gam = gam0 + 1;

        while abs(gam - gam0) > tol
            gam0 = gam;
            gam = ((1-mu)*(gam0+1)^2 / (1+2*mu+gam0*(2+mu+gam0)))^(1/3);
        end
        gamma(i) = gam;
        xL(i) = -mu - gam;

    end %if


    % solve for Jacobi constant of Lagrange point
    L(i) = energy([xL(i) 0 0 0],mu);
end %for
end % function
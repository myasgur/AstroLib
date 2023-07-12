function Cvals = cr3bp_contours(fig, mustar, localneck)
%Plot CR3BP zero-velocity curves corresponding to Lagrange points for given
%\mu^\star
%
%Inputs:
%   fig  (1)   Figure number or handle to plot in
%   mustar (1) \mu^\star value equal to m_2/(m_2 + m_1) where m_2 is the
%               less massive primary.
% localneck (1) Zoom in figure to neck region of ZVC around L1 (1), L2 (2)
%                  or view whole ZVC (0) ... or both L1 and L2 (12)
%
%Outputs:
%   Cvals  (4x1)  Jacobi constant values for L_1...4 in canonical units
%                 i.e., n = G = 1
Lpt = [1 2 3];
[~,~,xL123] = findGamma(Lpt,mustar); % x-positions of collinear Lagrange points
xL4         = [0.5-mustar, -sqrt(3)/2, 0, 0]; % full state of L4/5 points
x1234       = [xL123' zeros(length(xL123),3); xL4];

Cvals = energy(x1234,mustar);

if nargin < 3 || localneck==0
    tmpx = linspace(-1.5,1.5,1000);
    tmpy = tmpx;
elseif localneck == 1 || localneck == 2
    xL123 = xL123(localneck);
    offset = 2*10^(-ceil(-0.5*log10(mustar)));
    tmpx = linspace(xL123-offset,xL123+offset,1000);
    tmpy = linspace(-offset,offset,1000);
else
    xL1 = xL123(1); xL2 = xL123(2);
    span = (xL2-xL1)/2;
    tmpx = linspace(xL1-span/2,xL2+span/2,1000);
    tmpy = linspace(-span,span,1000);
end
[X,Y] = meshgrid(tmpx,tmpy);
r1 = sqrt((X+mustar).^2+Y.^2);
r2 = sqrt((X - (1-mustar)).^2+Y.^2);
U = -(X.^2+Y.^2)/2 - ((1-mustar)./r1 + mustar./r2 + 0.5*mustar*(1-mustar));

figure(fig)
hold on

if nargin<3 || localneck == 0
    contourf(X,Y,U,Cvals)
    plot([-mustar,1-mustar],[0,0],'r.','MarkerSize',15)
    plot(xL123,[0,0,0],'k.','MarkerSize',10)
    plot([xL4(1) xL4(1)],[xL4(2),-xL4(2)],'k.','MarkerSize',10)
else
    clevels = [Cvals(2),Cvals(2)].*(1-10^-ceil(-log10(mustar)));
    contourf(X,Y,U,clevels)
    if localneck == 1 || localneck == 2
        plot(xL123,0,'k.','MarkerSize',10)
    else
        plot(xL1,0,'k.',xL2,0,'k.',1-mustar,0,'r.','MarkerSize',10)
    end
end
colormap('winter')
axis('equal')
xlabel('$\mathbf{\hat{e}}_r$','Interpreter','LaTex')
ylabel('$\mathbf{\hat{e}}_\theta$','Interpreter','LaTex')
set(gca,'FontName','Times','FontSize',16)
function cr3bp_plot(t,X,mu,Lpt,fig)

n=size(X,2);
showEnergySurf = 12; % choose between 0, 1, 2, or 12
[L,~,xL] = findGamma(Lpt,mu);
 % construct state of Lagrange point
XL = [xL zeros(1,n-1)];

figure(fig)
tiledlayout(2,2,"Padding","compact")

switch n
    case 4
        % unpack state
        x = X(:,1); y = X(:,2); vx = X(:,3); vy = X(:,4);

        nexttile(1)
        plot(x,y,'k',x(1),y(1),'mo',XL(1),XL(2),'kx')
        % hold on
        % plot(-mu,0,'b*',1-mu,0,'k*')
        % hold off

        title('Rotating Frame');
        xlabel('x'); ylabel('y');
        axis equal; grid on; box on;

        % plot again in separate figure
        Cvals = cr3bp_contours(fig+1,mu,showEnergySurf);
        plot(x,y,'k',x(1),y(1),'mo',XL(1),XL(2),'kx')
        hold off
        set(gcf,"WindowStyle","docked")

        Xin = rot2iner(X,t,1,mu);
        [Xm1i,Xm2i] = rot2iner_prim(t,1,mu);

        figure(fig)
        nexttile(3)
        plot(Xin(:,1),Xin(:,2),'k',...
            Xin(1,1),Xin(1,2),'mo',...
            Xm1i(:,1),Xm1i(:,2),'b',...
            Xm2i(:,1),Xm2i(:,2),'r');
        title('Inertial Frame (m1-centered)');
        xlabel('x_{iner}'); ylabel('y_{iner}');
        axis equal; grid on; box on;

        Xin = rot2iner(X,t,2,mu);
        [Xm1i,Xm2i] = rot2iner_prim(t,2,mu);

        nexttile(4)
        plot(Xin(:,1),Xin(:,2),'k',...
            Xin(1,1),Xin(1,2),'mo',...
            Xm1i(:,1),Xm1i(:,2),'b',...
            Xm2i(:,1),Xm2i(:,2),'r');
        title('Inertial Frame (m2-centered)');
        xlabel('x_{iner}'); ylabel('y_{iner}');
        axis equal; grid on; box on;

    case 6
        % unpack state
        x = X(:,1); y = X(:,2); z = X(:,3); 
        vx = X(:,4); vy = X(:,5); vz = X(:,6);

        nexttile(1)
        plot3(x,y,z,'k',x(1),y(1),z(1),'mo',XL(1),XL(2),XL(3),'kx')
        title('Rotating Frame');
        xlabel('x'); ylabel('y');
        axis equal; grid on; box on;

        % plot again in separate figure with pseudopotential (planar)
        Cvals = cr3bp_contours(fig+1,mu,showEnergySurf);
        plot(x,y,'k',x(1),y(1),'mo',XL(1),XL(2),'kx')
        hold off
        set(gcf,"WindowStyle","docked")

        Xin = rot2iner(X,t,1,mu);
        [Xm1i,Xm2i] = rot2iner_prim(t,1,mu,0,1);

        figure(fig)
        nexttile(3)
        plot3(Xin(:,1),Xin(:,2),Xin(:,3),'k',...
              Xin(1,1),Xin(1,2),Xin(1,3),'mo',...
              Xm1i(:,1),Xm1i(:,2),Xm1i(:,3),'b',...
              Xm2i(:,1),Xm2i(:,2),Xm2i(:,3),'r');
        title('Inertial Frame (m1-centered)');
        xlabel('x_{iner}'); ylabel('y_{iner}'); zlabel('z_{iner}');
        axis equal; grid on; box on;

        Xin = rot2iner(X,t,2,mu);
        [Xm1i,Xm2i] = rot2iner_prim(t,2,mu,0,1);

        nexttile(4)
        plot3(Xin(:,1),Xin(:,2),Xin(:,3),'k',...
              Xin(1,1),Xin(1,2),Xin(1,3),'mo',...
              Xm1i(:,1),Xm1i(:,2),Xm1i(:,3),'b',...
              Xm2i(:,1),Xm2i(:,2),Xm2i(:,3),'r');
        title('Inertial Frame (m2-centered)');
        xlabel('x_{iner}'); ylabel('y_{iner}'); zlabel('z_{iner}');
        axis equal; grid on; box on;


end

% lastly, plot the energy to check the numerical stability
E = energy(X,mu);
nexttile(2)
hold on
plot(t,E,'k','LineWidth',0.75)
% plot([t(1) t(end)],[L L],'r--')
hold off
title('Jacobi Constant')
xlabel('t, time (2\pi = period of primaries)')
ylabel('E')
grid on; box on;

set(gcf,'WindowStyle','Docked')



end
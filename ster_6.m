
function ster_6
clear all; 
close all;
grid on;
tspan = [0 200]

x0 = [1.5 0.0 0.0 7.0];

% u = [0:0.1:50]
u=40
u_max = 0;
J_max = 0;
tf_max = 0;

% for i=1:length(u)        

        options_ode = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@eventFun);
        [tsol,xsol] = ode45(@(t,x) penicilin_production(t,x,u), tspan, x0, options_ode);

        x2 = xsol(:,2);
        x4 = xsol(:,4);

        J = x2(end)*x4(end)
%         J_tab(i) = J
%         if J > J_max
%             J_max = J;
%             u_max = u(i);
%         end
% end
% 
% J_max
% u_max
% plot(u,J_tab)
hold on;
plot(tsol,xsol)
%plot(14,J,'*')
title(['Trajektorie zmiennych stanu, u= ', num2str(u), ', J= ' , num2str(J)])

ylabel('Stê¿enie')

xlabel('Czas')

legend('x1(t)','x2(t)','x3(t)','x4(t)','J')
end

function [check,isterminal,direction] = eventFun(tsol,xsol)

check = ((xsol(1)& xsol(3)& xsol(4))>=0 & xsol(1)<=40 & xsol(3)<=25 & xsol(4)<=10);
isterminal = 1;
direction = [];

end

function dx = penicilin_production(t,x,u)

dx = zeros(4,1);
h1 = 0.11*(x(3)/(0.006*x(1)+x(3)));
h2 = 0.0055*(x(3)/(0.0001+x(3)*(1+10*x(3))));
h = [h1 h2];
dx(1) = h(1)*x(1)-u*(x(1)/(500*x(4)));
dx(2) = h(2)*x(1)-0.01*x(2)-u*(x(2)/(500*x(4)));
dx(3) = -h(1)*x(1)/0.47-h(2)*x(1)/1.2-x(1)*(0.029*x(3)/(0.0001+x(3)))+(u/x(4))*(1-x(3)/500);
dx(4) = u/500;

end
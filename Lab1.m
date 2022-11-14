function Lab1


tspan = [0.0  1.0];


x0 = [1.0  0.0];

 % założona wartość sterowania u(t) = [0  1]

u = 0.6;

 

options_ode = odeset('RelTol',1e-6,'AbsTol',1e-6);

[tsol,xsol] = ode45(@(t,x) catalyst_mixing_problem(t,x,u), tspan, x0, options_ode);

 

x3 = 1.0 - xsol(:,1) - xsol(:,2);

 

wskaznik_jakosci = x3(end)

 

 plot(tsol,xsol,tsol,x3)
title('Trajektorie zmiennych stanu')
ylabel('Stężenie')
xlabel('Czas')
legend('x1(t)','x2(t)','x3(t)')
grid on
grid minor

end


function dx = catalyst_mixing_problem(t,x,u)

dx = zeros(2,1);

dx(1) = u * (  10*x(2) - x(1)  );

dx(2) = u * (  x(1) - 10*x(2)  ) - ( 1 - u ) * x(2) ;

end
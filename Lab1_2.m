wskaznik_jakosci = zeros(1,100);
i = 1;
for u = 0:0.01:1
    
    tspan = [0.0  1.0];
    x0 = [1.0  0.0];
    options_ode = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [tsol,xsol] = ode45(@(t,x) catalyst_mixing_problem(t,x,u), tspan, x0, options_ode);
    x3 = 1.0 - xsol(:,1) - xsol(:,2);
    wskaznik_jakosci(i) = x3(end);
    i = i + 1;
end

u1 = 0:0.01:1;

plot(u1, wskaznik_jakosci);
title('Wykres zależności wartości wskaźnika jakości od wartości funkcji sterującej')
ylabel('Wskaźnik jakości')
xlabel('Funkcja sterująca')
grid on
grid minor

function dx = catalyst_mixing_problem(t,x,u)

dx = zeros(2,1);

dx(1) = u * (  10*x(2) - x(1)  );

dx(2) = u * (  x(1) - 10*x(2)  ) - ( 1 - u ) * x(2) ;

end
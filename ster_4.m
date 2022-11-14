function ster_4

 

% u = [ 0.9 0.9 0.9  0.088 0.088  0.0009  0.0009   ]

  u = [1.1 0.6 0.1  0.9 0.9  0.1  0.1];

wsk_jakosci = model_procesu(u)

 

end

 

 

function wskaznik_jakosci = model_procesu(u)

 

% model procesu po zastosowaniu metody strzałów wielopunktowych

% wybrano podział na 3 podprzedziały

 

tspan_1 = [0.0  1.0]/3;

tspan_2 = max(tspan_1) + [0.0  1.0]/3;

tspan_3 = max(tspan_2) + [0.0  1.0]/3;

 

x0_1 = [ 1.0   0.0  ];

x0_2 = [ u(4)  u(6) ];

x0_3 = [ u(5)  u(7) ];

 

options_ode = odeset('RelTol',1e-6,'AbsTol',1e-6);

 

[tsol_1,xsol_1] = ode45(@(t,x) catalyst_mixing_problem( t,x,u(1) ), tspan_1, x0_1, options_ode);

[tsol_2,xsol_2] = ode45(@(t,x) catalyst_mixing_problem( t,x,u(2) ), tspan_2, x0_2, options_ode);

[tsol_3,xsol_3] = ode45(@(t,x) catalyst_mixing_problem( t,x,u(3) ), tspan_3, x0_3, options_ode);

 

x3_1 = 1.0 - xsol_1(:,1) - xsol_1(:,2);

x3_2 = 1.0 - xsol_2(:,1) - xsol_2(:,2);

x3_3 = 1.0 - xsol_3(:,1) - xsol_3(:,2);

 

r1 = sum( (xsol_1(end, : ) - x0_2).^2 );

r2 = sum( (xsol_2(end, : ) - x0_3).^2 );

 

wskaznik_jakosci = x3_3(end) - 10^4*( r1 + r2 );

 

plot(tsol_1,xsol_1,tsol_1,x3_1, tsol_2, xsol_2,tsol_2,x3_2, tsol_3,xsol_3,tsol_3,x3_3 )

title('Trajektorie zmiennych stanu')

ylabel('Stężenie')

xlabel('Czas')

legend('x1_1(t)','x2_1(t)','x3_1(t)', 'x1_2(t)','x2_2(t)','x3_2(t)', 'x1_3(t)','x2_3(t)','x3_3(t)' )
%     figure
%     plot(tsol_1,xsol_1,tsol_1,x3_1)
%     title('Trajektorie zmiennych stanu')
%     ylabel('Stezenie')
%     xlabel('Czas')
%     legend('x1_1(t)','x2_1(t)','x3_1(t)', 'location', 'west');
%     figure
%     plot(tsol_2,xsol_2,tsol_2,x3_2)
%     title('Trajektorie zmiennych stanu')
%     ylabel('Stezenie')
%     xlabel('Czas')
%     legend('x1_2(t)','x2_2(t)','x3_2(t)', 'location', 'west');  
%     figure
%     plot(tsol_3,xsol_3,tsol_3,x3_3)
%     title('Trajektorie zmiennych stanu')
%     ylabel('Stezenie')
%     xlabel('Czas')
%     legend('x1_3(t)','x2_3(t)','x3_3(t)', 'location', 'west');
 

end

 

function dx = catalyst_mixing_problem(t,x,u)

dx = zeros(2,1);

dx(1) = u * (  10*x(2) - x(1)  );

dx(2) = u * (  x(1) - 10*x(2)  ) - ( 1 - u ) * x(2) ;

end
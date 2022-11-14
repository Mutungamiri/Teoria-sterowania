function ster_3_gotowe

 % zmienne decyzyjne

 tx = 170;    % wartoœæ z przedzia³u 1.0 - 179.0 min

T = 300;    % wartoœæ z przedzia³u 298 - 398 K

S = 0.05;   % wartoœæ z przedzia³u 0.0 - 0.1 m^3
maxwsk=0;
maxu=[0 0 0];
for tx =1:179
    for T = 298:398
        for S = 0:0.1:1
            u= [ tx T S ];



            % parametry uk³adu

            V1 = 0.1;      % m^3

            CA_1_0 = 2000; % mol/(m^3)

            CB0 = 600;     % mol/(m^3)



            % symulacja procesu w Reaktorze 1

            % zmienne stanu x = [CA CB CC CD CE CF]



            tspan_1 = [0.0 u(1)];

            x0_1 = [ CA_1_0  0.0  0.0  0.0  0.0  0.0 ];



            options_ode = odeset('RelTol',1e-6,'AbsTol',1e-6);

            [tsol_1,xsol_1] = ode45(@(t,x) Stage_1(t,x,u(2)), tspan_1, x0_1, options_ode);



            % Mieszanie

            V2 = V1 + u(3);

            CA_2_0 = V1*xsol_1(end,1) / V2 ;

            CB_2_0 = ( V1*xsol_1(end,2) + u(3)*CB0 ) / V2 ;

            CC_2_0 = V1*xsol_1(end,3) / V2 ;



            % symulacja procesu w Reaktorze 2

            tspan_2 = [u(1) 180];

            x0_2 = [ CA_2_0 CB_2_0 CC_2_0 0 0 0  ];

            [tsol_2,xsol_2] = ode45(@Stage_2, tspan_2, x0_2, options_ode);



            %wskaznik_jakosci

            wskaznik_jakosci = V2*xsol_2(end,4);
            if wskaznik_jakosci> maxwsk
                maxwsk = wskaznik_jakosci
                maxu=u
            end 


            subplot(1,2,1)

            %plot(tsol_1,xsol_1(:,1:3),tsol_2,xsol_2(:,1:3))

            title('Trajektorie zmiennych stanu CA(t), CB(t), CC(t)')

            ylabel('Stê¿enie')

            xlabel('Czas')

            legend('CA(t) - reaktor 1','CB(t) - reaktor 1','CC(t) - reaktor 1','CA(t) - reaktor 2','CB(t) - reaktor 2','CC(t) - reaktor 2')



            subplot(1,2,2)

            %plot(tsol_2,xsol_2(:,4:6))

            title('Trajektorie zmiennych stanu CD(t), CE(t), CF(t)')

            ylabel('Stê¿enie')

            xlabel('Czas')

            legend('CD(t)','CE(t)','CF(t)')
        end
    end
end
maxwsk
maxu
end

 

function dx = Stage_1(t,x,u)

CA = x(1); CB = x(2); CC = x(3);

CD = x(4); CE = x(5); CF = x(6);

Temp = u;

dx = zeros(6,1);

k1 = 0.0444 * exp( -2500/(Temp) );

k2 = 6889.0 * exp( -5000/(Temp) );

dx(1) = -2*k1*CA*CA;

dx(2) = k1*CA*CA - k2*CB;

dx(3) = k2*CB;

dx(4) = 0.0;

dx(5) = 0.0;

dx(6) = 0.0;

end

 

function dx = Stage_2(t,x)

CA = x(1); CB = x(2); CC = x(3);

CD = x(4); CE = x(5); CF = x(6);

dx = zeros(6,1);

dx(1) = 0.0;

dx(2) = -0.02*CB-0.05*CB-2*4.0*(10^(-5))*CB*CB;

dx(3) = 0.0;

dx(4) = 0.02*CB;

dx(5) = 0.05*CB;

dx(6) = 4.0*(10^(-5))*CB*CB;

end
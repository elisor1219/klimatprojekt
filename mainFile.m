%%Denna fil inh�ller l�sningar till klimatprojektet

%% Uppgift 1 
clc;clear;clf
%Svar p� fr�ga: Tror de skilljer sig �t, delvis f�r att de tv� modellerna
%anv�nder olika startv�rden, samt att beta kanske �r f�r�ndrad.
%Sedan verkar det som att RCP45 slutar vid runt �r 2150 medan v�r modell
%fors�tter till �r 2500.

run('utslappRCP45.m')

run('koncentrationerRCP45.m')

years = 1765:2500;
B = modellKolcykeln(years, 0.35, CO2Emissions);

subplot(2,1,1)
plot(years,B')
title('Our modell')
grid on

subplot(2,2,3)
plot(years,B(1,:)'*0.469)
title('CO_2 conc our modell')
grid on

subplot(2,2,4)
plot(years,CO2ConcRCP45)
title('CO_2 conc RCP45')
grid on

%% Uppgift 2
clc;clear;clf

run('utslappRCP45.m')

years = 1765:2500;
numberOfLines = 5;
beta = linspace(0.1,1,numberOfLines);
%Fixar en matris f�r koldioxidkoncentrationen f�r olika beta
B = zeros(3,length(years));
B_koldioxidkoncentrationen = zeros(numberOfLines, length(years));
B_biomassa = zeros(numberOfLines, length(years));
B_marken = zeros(numberOfLines, length(years));
for i=1:numberOfLines
    B = modellKolcykeln(years, beta(i), CO2Emissions);
    B_koldioxidkoncentrationen(i,:) = B(1,:)*0.469;
    B_biomassa(i,:) = B(2,:);
    B_marken(i,:) = B(3,:);
end

%Fixar s� att legend f�r r�tt beta-v�rde
legendPlot = string(beta);
for i=1:numberOfLines
    legendPlot(i) = append('\beta = ', legendPlot(i));
end

%Plottar koldioxidkoncentrationen
subplot(3,1,1)
plot(years,B_koldioxidkoncentrationen')
title('Koldioxidkoncentrationen')
grid on
legend(legendPlot,'Location','northwest')

%Plottar kol i biomassa
subplot(3,1,2)
plot(years,B_biomassa')
title('Kol i biomassa')
grid on
legend(legendPlot,'Location','northwest')

%Plottar kol i marken
subplot(3,1,3)
plot(years,B_marken')
title('Kol i marken')
grid on
legend(legendPlot,'Location','northwest')

%TODO KOlla fr�gan!!!
%------------------------------------------------------------------------------------------------------------------------------------------------------










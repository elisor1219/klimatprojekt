%%Denna fil inh?ller l?sningar till klimatprojektet

%% Uppgift 1 
clc;clear;clf
%Svar p? fr?ga: Tror de skilljer sig ?t, delvis f?r att de tv? modellerna
%anv?nder olika startv?rden, samt att beta kanske ?r f?r?ndrad.
%Sedan verkar det som att RCP45 slutar vid runt ?r 2150 medan v?r modell
%fors?tter till ?r 2500.

run('koncentrationerRCP45.m')

Year = 1765:2500;
B = modellKolcykeln(Year, 0.35);

fontSize = 15;

subplot(2,1,1)
plot(Year,B')
title('Our modell', 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('GtC (Gigaton kol)', 'FontSize', fontSize)
legend(["Atmosf?r";"Ovan mark";"Under mark"],'Location','northwest', 'FontSize', fontSize)
grid on

subplot(2,2,3)
plot(Year,B(1,:)'*0.469)
title('CO_2 conc our modell', 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('ppm', 'FontSize', fontSize)
grid on

subplot(2,2,4)
plot(Year,CO2ConcRCP45)
title('CO_2 conc RCP45', 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('ppm', 'FontSize', fontSize)
grid on

%% Uppgift 2
clc;clear;clf
%Svar p? fr?ga: %H?gre beta v?rde ger mer kol i marken och biomassan,
%d?remot sjunker koldioxidkoncentrationen i atmosf?ren. 
%Detta ?r d? vi ?kar beta, vilket ?r CO_2-fertiliseringsfaktorn, som 
%kommer ?ka ?verf?rningen av CO_2 till biomassan som i sin tur ?kar 
%?verf?rningen till marken.


Year = 1765:2500;
numberOfLines = 5;
beta = linspace(0.1,1,numberOfLines);
%Fixar en matris f?r koldioxidkoncentrationen f?r olika beta
B_koldioxidkoncentrationen = zeros(numberOfLines, length(Year));
B_biomassa = zeros(numberOfLines, length(Year));
B_marken = zeros(numberOfLines, length(Year));
for i=1:numberOfLines
    B = modellKolcykeln(Year, beta(i));
    B_koldioxidkoncentrationen(i,:) = B(1,:)*0.469;
    B_biomassa(i,:) = B(2,:);
    B_marken(i,:) = B(3,:);
end

%Fixar s? att legend f?r r?tt beta-v?rde
legendPlot = string(beta);
for i=1:numberOfLines
    legendPlot(i) = append('\beta = ', legendPlot(i));
end

fontSize = 15;

%Plottar koldioxidkoncentrationen
subplot(3,1,1)
plot(Year,B_koldioxidkoncentrationen')
title('Koldioxidkoncentrationen', 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('ppm', 'FontSize', fontSize)
grid on
legend(legendPlot,'Location','northwest')

%Plottar kol i biomassahttps://github.com/elisor1219/klimatprojekt
subplot(3,1,2)
plot(Year,B_biomassa')
title('Kol i biomassa', 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('ppm', 'FontSize', fontSize)
grid on
legend(legendPlot,'Location','northwest')

%Plottar kol i marken
subplot(3,1,3)
plot(Year,B_marken')
title('Kol i marken', 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('ppm', 'FontSize', fontSize)
grid on
legend(legendPlot,'Location','northwest')

%% Uppgift 3
clc;clear;clf

U = [0 140 560 1680];

%Plottar
Year = 1:500;
plotTemp = zeros(length(Year),4);
for i = 1:length(Year)
    plotTemp(i,1) = impulssvaret(i,1,U(1));
    plotTemp(i,2) = impulssvaret(i,1,U(2));
    plotTemp(i,3) = impulssvaret(i,1,U(3));
    plotTemp(i,4) = impulssvaret(i,1,U(4));
end

%Fixar s? att legend f?r r?tt U-v?rde
legendPlot = string(U);
for i=1:4
    legendPlot(i) = append(legendPlot(i), ' GtC');
end

fontSize = 15;


plot(Year,plotTemp)
title('Impulsrespons f?r CO_2 beroende p? tidigare kumulativa utsl?pp', 'FontSize', fontSize)
xlabel('?r efter utsl?ppspuls', 'FontSize', fontSize)
ylabel('Andel kvar i atmosf?ren', 'FontSize', fontSize)
grid on
legend(legendPlot,'Location','northeast', 'FontSize', fontSize)

%% Uppgift 4
clc;clear;clf

run('koncentrationerRCP45.m')

%St?ller upp utsl?ppet
run('utslappRCP45.m')
U = CO2Emissions;

Year = 1765:2500;
M = tidsdeskretFaltning(Year, U);



%Plottar
fontSize = 14;


plot(Year, M*0.469)
title('Our new modell', 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('ppm', 'FontSize', fontSize)
grid on
hold on
plot(Year,CO2ConcRCP45)
hold off
legendPlot = ["CO_2 conc our modell"
              "CO_2 conc RCP45"];
legend(legendPlot,'Location','northwest', 'FontSize', fontSize)

%% Uppgift 5
clf

fullFileName = '../figures/boxmodell.png';
if ~isfile(fullFileName)
  error('Warning: file does not exist and is not necessary, just go to the next assignment:\n%s',fullFileName);
end

boxmodell = importdata(fullFileName);
image(boxmodell)

%% Uppgift 6
clc;clear;clf

run('koncentrationerRCP45.m')

Year = 1765:2500;
beta = 0.28;
newModell = modellSammankopplade(Year, beta);

%Plottar
fontSize = 15;

subplot(2,1,1)
plot(Year,newModell')
titleString = append('Two models in one (\beta = ', string(beta),')');
title(titleString, 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('GtC (Gigaton kol)', 'FontSize', fontSize)
legend(["Atmosf?r";"Ovan mark";"Under mark";"Hav"],'Location','west', 'FontSize', 11)
grid on

subplot(2,1,2)
plot(Year, newModell(1,:)*0.469)
titleString = append('Koncentration av CO_2 i atmosf?ren (\beta = ', string(beta),')');
title(titleString, 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('ppm', 'FontSize', fontSize)
grid on
hold on
plot(Year,CO2ConcRCP45(1:length(Year)))
hold off
legendPlot = ["CO_2 conc our modell"
              "CO_2 conc RCP45"];
legend(legendPlot,'Location','southeast', 'FontSize', fontSize)

%% Uppgift 7
%?ver sikt f?rflyttar sig mycket CO_2 till haven.
%Mellan 1765-2100 ?r det atmosf?ren som ?ndras mest, sedan f?ljder
%ovan mark, under mark och haven ganska lika kurvor.

%-%-%-Beta = 0.28*2 & k = 3.06E-3-%-%-%
%Av att dubbla CO_2-fertiliseringsfaktorn s? ?kas antal CO_2 som f?rflyttas
%till biomassa ovan och under mark. M?ngden CO_2 som f?rflyttar sig till
%havet f?rblir hyfsat of?r?ndrat.

%-%-%-Beta = 0.28 & k = 3.06*2E-3-%-%-%
%Dubblar vi k ?kar vi hur snabbt 
%havet blir m?ttat s? mindre kol kommer f?rlytta sig till havet och
%mer kol kommer befinna sig i atmosf?ren samt biomassa ovan och under mark

%-%-%-Beta = 0.28/2 & k = 3.06E-3-%-%-%
%Om vi halverar beta minskas hur mycket kol som kommer f?rflytta sig 
%till biomassa ovan och under mark. Detta g?r att mer kol kommer samla
%sig i atmosf?ren. M?ngden CO_2 som f?rflyttar sig till havet f?rblir 
%hyfsat of?r?ndrat.

%-%-%-Beta = 0.28 & k = 3.06/2E-3-%-%-%
%Halverar vi k ?kar vi hur mycket CO_2 havet kan ta upp innan det blir
%m?ttat. Detta g?r att mer CO_2 kommer f?rlytta sig fr?n atmosf?ren
%till havet under kortare tid. 

%% Uppgift 8
clc;clear;clf

run('koncentrationerRCP45.m');
run('radiativeForcingRCP45.m');


Year = 1765:2500;
RF = RF_CO2(CO2ConcRCP45);

fontSize = 15;
plot(Year, RF)
titleString = append('Radiative forcing f?r CO_2');
title(titleString, 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('W/M^2', 'FontSize', fontSize)
grid on
hold on
plot(Year,CO2RadForc)
hold off
legendPlot = ["Ber?knad RF-C0_2"
              "Data fr?n RCP45"];
legend(legendPlot,'Location','southeast', 'FontSize', fontSize)

%% Uppgift 9
clc;clear;clf

run('radiativeForcingRCP45.m');

s = 1;
totalRadForceExclCO2 = totRadForcExclCO2AndAerosols + totRadForcAerosols * s;

Year = 1765:2500;

fontSize = 15;
plot(Year, totalRadForceExclCO2)
titleString = append('Radiative forcing utan CO_2');
title(titleString, 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('W/M^2', 'FontSize', fontSize)
grid on

%% Uppgift 10a
%Svar p? fr?ga: Runt ?r 3321 s? ?r deltaT_1 = deltaT_2 = RF*lamda.
%R?knat p? en tolerans p? 0.005.
clc;clear;clf

time = 1:10000;
RF = ones(length(time),1);
RF(1) = 0; %Enligt uppgiften
lambda = 0.9; %spann p? 0.5-1.3
k = 0.6; %spann p? 0.2-1

deltaT = vattenModell(time, RF, lambda, k);
plot(time, deltaT(:,:))

tol = 0.005;
for i = 2:length(time)
   if abs(deltaT(1,i) - deltaT(2,i)) <= tol && abs(deltaT(1,i) - RF(i)*lambda) <= tol
       disp(['?r ', num2str(i), ' ?r deltaT_1 = deltaT_2 = RF*lamda'])
       break
   end
    
end 

%% Uppgift 10b
clc;clear;clf
%-%-%-Lamda = 0.5 & k = 0.6-%-%-%
%?r 6 n?r deltaT_1 e-folding time
%?r 592 n?r deltaT_2 e-folding time
%?r 4024 ?r deltaT_1 = deltaT_2 = RF*lamda

%-%-%-Lamda = 1.3 & k = 0.6-%-%-%
%?r 154 n?r deltaT_1 e-folding time
%?r 814 n?r deltaT_2 e-folding time
%?r 2553 ?r deltaT_1 = deltaT_2 = RF*lamda

%-%-%-Lamda = 0.9 & k = 0.2-%-%-%
%?r 9 n?r deltaT_1 e-folding time
%?r 1606 n?r deltaT_2 e-folding time
%?r 8041 ?r deltaT_1 = deltaT_2 = RF*lamda

%-%-%-Lamda = 0.9 & k = 0.1-%-%-%
%?r 140 n?r deltaT_1 e-folding time
%?r 523 n?r deltaT_2 e-folding time
%?r 2355 ?r deltaT_1 = deltaT_2 = RF*lamda

time = 1:10000;
RF = ones(length(time),1);
RF(1) = 0; %Enligt uppgiften
lambda = 0.9; %spann p? 0.5-1.3
k = 0.6; %spann p? 0.2-1

deltaT = vattenModell(time, RF, lambda, k);
plot(time, deltaT(:,:))

tol = 0.005;
e_foldingTime = (1-exp(-1))*(RF(end)*lambda);
foldingTimeDone = [false, false];
for i = 2:length(time)
   if abs(deltaT(1,i) - deltaT(2,i)) <= tol && abs(deltaT(1,i) - RF(i)*lambda) <= tol
       disp(['?r ', num2str(i), ' ?r deltaT_1 = deltaT_2 = RF*lamda'])
       break
   end
   
   if deltaT(1,i) >= e_foldingTime && ~foldingTimeDone(1)
       foldingTimeDone(1) = true;
       disp(['?r ', num2str(i), ' n?r deltaT_1 e-folding time'])
   end
   
   if deltaT(2,i) >= e_foldingTime && ~foldingTimeDone(2)
       foldingTimeDone(2) = true;
       disp(['?r ', num2str(i), ' n?r deltaT_2 e-folding time'])
   end
end 

%% Uppgift 10c
clc;clear;clf
%Svar p? fr?ga: Japp, lagen om konservering av energi g?ller.
%L?gre lamda = mer enrgi ut i rymnden och mindre energi i systemet.
%H?gre lamda = tv?rtom fr?n ovan.
%L?gre kappa = mindre energi i haven = mer energi i rymden (konservering av energi).
%H?gre kappa = tv?rtom fr?n ovan.

time = 1:200;
RF = ones(length(time),1);
RF(1) = 0; %Enligt uppgiften
lambda = 0.9; %spann p? 0.5-1.3
k = 0.6; %spann p? 0.2-1

[deltaT, C_1, ~] = vattenModell(time, RF, lambda, k);

tol = 0.005;
e_foldingTime = (1-exp(-1))*(RF(end)*lambda);
foldingTimeDone = [false, false];
for i = 2:length(time)
   if abs(deltaT(1,i) - deltaT(2,i)) <= tol && abs(deltaT(1,i) - RF(i)*lambda) <= tol
       disp(['?r ', num2str(i), ' ?r deltaT_1 = deltaT_2 = RF*lamda'])
       break
   end
   
   if deltaT(1,i) >= e_foldingTime && ~foldingTimeDone(1)
       foldingTimeDone(1) = true;
       disp(['?r ', num2str(i), ' n?r deltaT_1 e-folding time'])
   end
   
   if deltaT(2,i) >= e_foldingTime && ~foldingTimeDone(2)
       foldingTimeDone(2) = true;
       disp(['?r ', num2str(i), ' n?r deltaT_2 e-folding time'])
   end
end 

%?r det energibevarande?
%Ackumulerad energi (i haven)
ackumuleradEnergiHaven = [sum(gradient(deltaT(1,:))*C_1,2),sum(k*(deltaT(1,:)-deltaT(2,:)),2)];
energiIn = sum(RF);
energiUt = sum(deltaT(1,:)/lambda);
disp(append('ackumuleradEnergiHaven = ', string(sum(ackumuleradEnergiHaven))))
disp(append('energiIn - energiUt = ', string(energiIn-energiUt)))


fontSize = 13;

plot(time, gradient(deltaT(1,:))*C_1)
hold on
plot(time, k*(deltaT(1,:)-deltaT(2,:)))
%plot(time, deltaT(2,:)*C_2)
plot(time, deltaT(1,:)/lambda)
plot(time, RF)
titleString = append('Energifl?det vid \lambda = ', string(lambda), ' och \kappa = ', string(k));
title(titleString, 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('W/M^2', 'FontSize', fontSize)
grid on
axis([0 200 0 1.05 ])
legendPlot = ["Ytbox energiupptag"
              "Djupbox energiupptag"
              "Energi ut i rymden"
              "Radiative forcing"];
legend(legendPlot,'Location','east', 'FontSize', fontSize)
hold off


%% Uppgift 11 abc
clc;clear;clf
%a) Svar p? fr?ga: Referencperioden p?verkar svaret med att f?rflytta
%                  modellen l?ngre upp eller ner.
%   Svar p? fr?ga: Referensperiod mellan 1766-1786 ?r bra f?r att se
%                  medeltemperatur?kning f?r f?rindustriell tid.
%b) Svar p? fr?ga: N?r lamda ?r stor ?kar temeraturen snabbare, tv?rtom
%                  annars.
%   Svar p? fr?ga: N?r k ?r stor ?r det mindre fluktueringar mellan st?rsta
%                  och minsta temperaturskillnad, tv?rtom annars.
%   Svar p? f?rga: lambda = 0.7, k = 0.8, s = 0.9
%                  lambda = 1.2, k = 1, s = 1.5
%                  lambda = 0.82, k = 0.89, s = 1.27
%   Svar p? fr?ga: Avg?ra n?r det ?r som minst os?kert med v?nderna och
%                  utg? fr?n det.

run('NASA_GISS.m')

Year = 1765:2019;
referensperiod = 1951:1980;
NASAYear = 1880:2019;

%Variabler
beta = 0.35; %spann 0.1-0.8
lambda = 0.82; %spann p? 0.5-1.3
k = 0.89; %spann p? 0.2-1
s = 1.27;   %spann olkart

[~, temperaturCase] = modellTemp(Year, beta, lambda, k, s);


globalMedeltemperaturCase = sum(temperaturCase,1);
globalMedeltemperaturKorigerd = globalMedeltemperaturCase - mean(globalMedeltemperaturCase(referensperiod-Year(1)+1));


diff = sum(abs(globalMedeltemperaturKorigerd(NASAYear - Year(1)+1) - TAnomali));
disp(append('Skillnad mellan NASA och v?r modell = ', string(diff)))


fontSize = 15;

plot(Year, globalMedeltemperaturKorigerd)
hold on
%NASA
plot(NASAYear,TAnomali)
titleString = append('Medeltemperatur?kning baserat p? medelv?rdet ?ver ', string(referensperiod(1)), ' - ', string(referensperiod(end)));
subtiTitle = append('(\beta = ', string(beta),', \lambda = ', string(lambda), ', \kappa = ', string(k),', s = ',string(s),')');
title([titleString, subtiTitle], 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('Medeltemperatur?kning', 'FontSize', fontSize)
legend(["Ber?knad medeltemperatur";"NASA medeltemperatur"],'Location','northwest', 'FontSize', fontSize)
hold off
axis([Year(1) Year(end) -2 2])
grid on

%% Uppgift 12 ab
clc;clf;clear

run('utslappRCP45.m')

Year = 1765:2100;
referensperiod = 1765:1865;

globalMedeltemperaturKorigerd = zeros(3,length(Year));

beta = 0.35; %spann 0.1-0.8
lambda = 0.8; %spann p? 0.5-1.3
k = 0.89; %spann p? 0.2-1
s = 1.27;   %spann olkart

%Case 1 -%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

[~, temperaturCase] = case1(Year, beta, lambda, k, s);


globalMedeltemperaturCase = sum(temperaturCase,1);
globalMedeltemperaturKorigerd(1,:) = globalMedeltemperaturCase - mean(globalMedeltemperaturCase(referensperiod-Year(1)+1));

%Case 2 -%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

[~, temperaturCase2] = case2(Year, beta, lambda, k, s);


globalMedeltemperaturCase2 = sum(temperaturCase2,1);
globalMedeltemperaturKorigerd(2,:) = globalMedeltemperaturCase2 - mean(globalMedeltemperaturCase2(referensperiod-Year(1)+1));

%Case 3 -%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

[~, temperaturCase3] = case3(Year, beta, lambda, k, s);


globalMedeltemperaturCase3 = sum(temperaturCase3,1);
globalMedeltemperaturKorigerd(3,:) = globalMedeltemperaturCase3 - mean(globalMedeltemperaturCase3(referensperiod-Year(1)+1));



fontSize = 15;

plot(Year, globalMedeltemperaturKorigerd)
hold on
titleString = append('Medeltemperatur?kning baserat p? medelv?rdet ?ver ', string(referensperiod(1)), ' - ', string(referensperiod(end)));
subtiTitle = append('(\beta = ', string(beta),', \lambda = ', string(lambda), ', \kappa = ', string(k),', s = ',string(s),')');
title([titleString, subtiTitle], 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('Temperatur', 'FontSize', fontSize)
legendPlot = ["Alternativ i"
              "Alternativ ii"
              "Alternativ iii"];
legend(legendPlot,'Location','northwest', 'FontSize', fontSize)
hold off
axis([Year(1) Year(end) -2 5])
grid on

%% Uppgift 12 c
clc;clf;clear

run('utslappRCP45.m')

Year = 1765:2100;
referensperiod = 1765:1865;

globalMedeltemperaturKorigerd = zeros(6,length(Year));

beta = 0.35; %spann 0.1-0.8
%lambda = 0.5; %spann p? 0.5-1.3
k = 0.89; %spann p? 0.2-1
s = 1.27;   %spann olkart

%Case 1 -%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

lambda = 0.5;
[~, temperaturCase] = case1(Year, beta, lambda, k, s);

globalMedeltemperaturCase = sum(temperaturCase,1);
globalMedeltemperaturKorigerd(1,:) = globalMedeltemperaturCase - mean(globalMedeltemperaturCase(referensperiod-Year(1)+1));

lambda = 1.3;
[~, temperaturCase] = case1(Year, beta, lambda, k, s);

globalMedeltemperaturCase = sum(temperaturCase,1);
globalMedeltemperaturKorigerd(4,:) = globalMedeltemperaturCase - mean(globalMedeltemperaturCase(referensperiod-Year(1)+1));

%Case 2 -%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

lambda = 0.5;
[~, temperaturCase] = case2(Year, beta, lambda, k, s);

globalMedeltemperaturCase = sum(temperaturCase,1);
globalMedeltemperaturKorigerd(2,:) = globalMedeltemperaturCase - mean(globalMedeltemperaturCase(referensperiod-Year(1)+1));

lambda = 1.3;
[~, temperaturCase] = case2(Year, beta, lambda, k, s);

globalMedeltemperaturCase = sum(temperaturCase,1);
globalMedeltemperaturKorigerd(5,:) = globalMedeltemperaturCase - mean(globalMedeltemperaturCase(referensperiod-Year(1)+1));

%Case 3 -%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

lambda = 0.5;
[~, temperaturCase] = case3(Year, beta, lambda, k, s);

globalMedeltemperaturCase = sum(temperaturCase,1);
globalMedeltemperaturKorigerd(3,:) = globalMedeltemperaturCase - mean(globalMedeltemperaturCase(referensperiod-Year(1)+1));

lambda = 1.3;
[~, temperaturCase] = case3(Year, beta, lambda, k, s);

globalMedeltemperaturCase = sum(temperaturCase,1);
globalMedeltemperaturKorigerd(6,:) = globalMedeltemperaturCase - mean(globalMedeltemperaturCase(referensperiod-Year(1)+1));




fontSize = 15;

%plot(Year, globalMedeltemperaturKorigerd(1:3,:))
hold on
%patch([Year,flip(Year)],[globalMedeltemperaturKorigerd(2,:),flip(globalMedeltemperaturKorigerd(5,:))], 'r', 'FaceAlpha', 0.3)
plot(Year,globalMedeltemperaturKorigerd(1:3,:))
plot(Year,globalMedeltemperaturKorigerd(4,:), '--', 'color', '[0, 0.4470, 0.7410]')
plot(Year,globalMedeltemperaturKorigerd(5,:), '--', 'color', '[0.8500, 0.3250, 0.0980]')
plot(Year,globalMedeltemperaturKorigerd(6,:), '--', 'color', '[0.9290, 0.6940, 0.1250]')
titleString = append('Medeltemperatur?kning baserat p? medelv?rdet ?ver ', string(referensperiod(1)), ' - ', string(referensperiod(end)));
subtiTitle = append('(\beta = ', string(beta),', \lambda = 0.5 - 1.3 \kappa = ', string(k),', s = ',string(s),')');
title([titleString, subtiTitle], 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('Temperatur', 'FontSize', fontSize)
legendPlot = ["Alternativ i   low"
              "Alternativ ii  low"
              "Alternativ iii low"
              "Alternativ i   high"
              "Alternativ ii  high"
              "Alternativ iii high"];
legend(legendPlot,'Location','northwest', 'FontSize', fontSize)
hold off
axis([Year(1) Year(end) -2 5])
grid on

%% Uppgift 12d se text
%% Uppgift 13a 

clc;clf;clear

run('utslappRCP45.m')

Year = 1765:2200;
referensperiod = 1765:1865;

beta = 0.35; %spann 0.1-0.8
lambda = 1.3; %spann p? 0.5-1.3
k = 0.89; %spann p? 0.2-1
s = 1.27;   %spann olkart

RF_geoEngineering = zeros(1,length(Year));
%S?tter 2050 - 2100 att vi minskar med 4 W/m^2
RF_geoEngineering(286:336) = 4;

[deltaT, temperaturCase, C_1] = case3(Year, beta, lambda, k, s, RF_geoEngineering);

globalMedeltemperaturCase = sum(temperaturCase,1);
globalMedeltemperaturKorigerd(1,:) = globalMedeltemperaturCase - mean(globalMedeltemperaturCase(referensperiod-Year(1)+1));

[~, temperaturCase] = case3(Year, beta, lambda, k, s);

globalMedeltemperaturCase = sum(temperaturCase,1);
globalMedeltemperaturKorigerd(2,:) = globalMedeltemperaturCase - mean(globalMedeltemperaturCase(referensperiod-Year(1)+1));




fontSize = 15;

hold on
fill([2050, 2100, 2100, 2050],[-10, -10, 20, 20],'g','FaceAlpha', 0.25, 'LineStyle','none')
plot(Year,globalMedeltemperaturKorigerd(2,:), 'color' ,'[0.6350, 0.0780, 0.1840]')
plot(Year,globalMedeltemperaturKorigerd(1,:), 'color' ,'[0.9290, 0.6940, 0.1250]')
titleString = append('Medeltemperatur?kning baserat p? medelv?rdet ?ver ', string(referensperiod(1)), ' - ', string(referensperiod(end)));
subtiTitle = append('(\beta = ', string(beta),', \lambda = ', string(lambda), ', \kappa = ', string(k),', s = ',string(s),')');
title([titleString, subtiTitle], 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('Temperatur', 'FontSize', fontSize)
legendPlot = ["?kar aerosoler i atmosf?ren"
              "Alternativ iii utan ingrepp"
              "Alternativ iii med extra aerosoler"];
legend(legendPlot,'Location','northwest', 'FontSize', fontSize)
hold off
axis([Year(1) Year(end) -2 10])
grid on

%% Uppgift 13b se text



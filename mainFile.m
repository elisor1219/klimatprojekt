%%Denna fil inhåller lösningar till klimatprojektet

%% Uppgift 1 
clc;clear;clf
%Svar på fråga: Tror de skilljer sig åt, delvis för att de två modellerna
%använder olika startvärden, samt att beta kanske är förändrad.
%Sedan verkar det som att RCP45 slutar vid runt år 2150 medan vår modell
%forsätter till år 2500.

run('koncentrationerRCP45.m')

Year = 1765:2500;
B = modellKolcykeln(Year, 0.35);

fontSize = 15;

subplot(2,1,1)
plot(Year,B')
title('Our modell', 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('GtC (Gigaton kol)', 'FontSize', fontSize)
legend(["Atmosfär";"Ovan mark";"Under mark"],'Location','northwest', 'FontSize', fontSize)
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
%Svar på fråga: %Högre beta värde ger mer kol i marken och biomassan,
%däremot sjunker koldioxidkoncentrationen i atmosfären. 
%Detta är då vi ökar beta, vilket är CO_2-fertiliseringsfaktorn, som 
%kommer öka överförningen av CO_2 till biomassan som i sin tur ökar 
%överförningen till marken.


Year = 1765:2500;
numberOfLines = 5;
beta = linspace(0.1,1,numberOfLines);
%Fixar en matris för koldioxidkoncentrationen för olika beta
B_koldioxidkoncentrationen = zeros(numberOfLines, length(Year));
B_biomassa = zeros(numberOfLines, length(Year));
B_marken = zeros(numberOfLines, length(Year));
for i=1:numberOfLines
    B = modellKolcykeln(Year, beta(i));
    B_koldioxidkoncentrationen(i,:) = B(1,:)*0.469;
    B_biomassa(i,:) = B(2,:);
    B_marken(i,:) = B(3,:);
end

%Fixar så att legend får rätt beta-värde
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

%Plottar kol i biomassa
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

%Fixar så att legend får rätt U-värde
legendPlot = string(U);
for i=1:4
    legendPlot(i) = append(legendPlot(i), ' GtC');
end

fontSize = 15;


plot(Year,plotTemp)
title('Impulsrespons för CO_2 beroende på tidigare kumulativa utsläpp', 'FontSize', fontSize)
xlabel('År efter utsläppspuls', 'FontSize', fontSize)
ylabel('Andel kvar i atmosfären', 'FontSize', fontSize)
grid on
legend(legendPlot,'Location','northeast', 'FontSize', fontSize)

%% Uppgift 4
clc;clear;clf

run('koncentrationerRCP45.m')

%Ställer upp utsläppet
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
legend(["Atmosfär";"Ovan mark";"Under mark";"Hav"],'Location','west', 'FontSize', 11)
grid on

subplot(2,1,2)
plot(Year, newModell(1,:)*0.469)
titleString = append('Koncentration av CO_2 i atmosfären (\beta = ', string(beta),')');
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
error("Nothing to run here, see source code for answers.")
%Över sikt förflyttar sig mycket CO_2 till haven.
%Mellan 1765-2100 är det atmosfären som ändras mest, sedan följder
%ovan mark, under mark och haven ganska lika kurvor.

%-%-%-Beta = 0.28*2 & k = 3.06E-3-%-%-%
%Av att dubbla CO_2-fertiliseringsfaktorn så ökas antal CO_2 som förflyttas
%till biomassa ovan och under mark. Mängden CO_2 som förflyttar sig till
%havet förblir hyfsat oförändrat.

%-%-%-Beta = 0.28 & k = 3.06*2E-3-%-%-%
%Dubblar vi k ökar vi hur snabbt 
%havet blir mättat så mindre kol kommer förlytta sig till havet och
%mer kol kommer befinna sig i atmosfären samt biomassa ovan och under mark

%-%-%-Beta = 0.28/2 & k = 3.06E-3-%-%-%
%Om vi halverar beta minskas hur mycket kol som kommer förflytta sig 
%till biomassa ovan och under mark. Detta gör att mer kol kommer samla
%sig i atmosfären. Mängden CO_2 som förflyttar sig till havet förblir 
%hyfsat oförändrat.

%-%-%-Beta = 0.28 & k = 3.06/2E-3-%-%-%
%Halverar vi k ökar vi hur mycket CO_2 havet kan ta upp innan det blir
%mättat. Detta gör att mer CO_2 kommer förlytta sig från atmosfären
%till havet under kortare tid. 







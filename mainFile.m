%%Denna fil inh�ller l�sningar till klimatprojektet

%% Uppgift 1 
clc;clear;clf
%Svar p� fr�ga: Tror de skilljer sig �t, delvis f�r att de tv� modellerna
%anv�nder olika startv�rden, samt att beta kanske �r f�r�ndrad.
%Sedan verkar det som att RCP45 slutar vid runt �r 2150 medan v�r modell
%fors�tter till �r 2500.

run('koncentrationerRCP45.m')

Year = 1765:2500;
B = modellKolcykeln(Year, 0.35);

fontSize = 15;

subplot(2,1,1)
plot(Year,B')
title('Our modell', 'FontSize', fontSize)
xlabel('Year', 'FontSize', fontSize)
ylabel('GtC (Gigaton kol)', 'FontSize', fontSize)
legend(["Atmosf�r";"Ovan mark";"Under mark"],'Location','northwest', 'FontSize', fontSize)
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
%Svar p� fr�ga: %H�gre beta v�rde ger mer kol i marken och biomassan,
%d�remot sjunker koldioxidkoncentrationen i atmosf�ren. 
%Detta �r d� vi �kar beta, vilket �r CO_2-fertiliseringsfaktorn, som 
%kommer �ka �verf�rningen av CO_2 till biomassan som i sin tur �kar 
%�verf�rningen till marken.


Year = 1765:2500;
numberOfLines = 5;
beta = linspace(0.1,1,numberOfLines);
%Fixar en matris f�r koldioxidkoncentrationen f�r olika beta
B_koldioxidkoncentrationen = zeros(numberOfLines, length(Year));
B_biomassa = zeros(numberOfLines, length(Year));
B_marken = zeros(numberOfLines, length(Year));
for i=1:numberOfLines
    B = modellKolcykeln(Year, beta(i));
    B_koldioxidkoncentrationen(i,:) = B(1,:)*0.469;
    B_biomassa(i,:) = B(2,:);
    B_marken(i,:) = B(3,:);
end

%Fixar s� att legend f�r r�tt beta-v�rde
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

%Fixar s� att legend f�r r�tt U-v�rde
legendPlot = string(U);
for i=1:4
    legendPlot(i) = append(legendPlot(i), ' GtC');
end

fontSize = 15;


plot(Year,plotTemp)
title('Impulsrespons f�r CO_2 beroende p� tidigare kumulativa utsl�pp', 'FontSize', fontSize)
xlabel('�r efter utsl�ppspuls', 'FontSize', fontSize)
ylabel('Andel kvar i atmosf�ren', 'FontSize', fontSize)
grid on
legend(legendPlot,'Location','northeast', 'FontSize', fontSize)

%% Uppgift 4
clc;clear;clf

run('koncentrationerRCP45.m')

%St�ller upp utsl�ppet
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
legend(["Atmosf�r";"Ovan mark";"Under mark";"Hav"],'Location','west', 'FontSize', 11)
grid on

subplot(2,1,2)
plot(Year, newModell(1,:)*0.469)
titleString = append('Koncentration av CO_2 i atmosf�ren (\beta = ', string(beta),')');
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
%�ver sikt f�rflyttar sig mycket CO_2 till haven.
%Mellan 1765-2100 �r det atmosf�ren som �ndras mest, sedan f�ljder
%ovan mark, under mark och haven ganska lika kurvor.

%-%-%-Beta = 0.28*2 & k = 3.06E-3-%-%-%
%Av att dubbla CO_2-fertiliseringsfaktorn s� �kas antal CO_2 som f�rflyttas
%till biomassa ovan och under mark. M�ngden CO_2 som f�rflyttar sig till
%havet f�rblir hyfsat of�r�ndrat.

%-%-%-Beta = 0.28 & k = 3.06*2E-3-%-%-%
%Dubblar vi k �kar vi hur snabbt 
%havet blir m�ttat s� mindre kol kommer f�rlytta sig till havet och
%mer kol kommer befinna sig i atmosf�ren samt biomassa ovan och under mark

%-%-%-Beta = 0.28/2 & k = 3.06E-3-%-%-%
%Om vi halverar beta minskas hur mycket kol som kommer f�rflytta sig 
%till biomassa ovan och under mark. Detta g�r att mer kol kommer samla
%sig i atmosf�ren. M�ngden CO_2 som f�rflyttar sig till havet f�rblir 
%hyfsat of�r�ndrat.

%-%-%-Beta = 0.28 & k = 3.06/2E-3-%-%-%
%Halverar vi k �kar vi hur mycket CO_2 havet kan ta upp innan det blir
%m�ttat. Detta g�r att mer CO_2 kommer f�rlytta sig fr�n atmosf�ren
%till havet under kortare tid. 







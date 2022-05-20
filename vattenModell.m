function [deltaT, C_1, C_2] = vattenModell(years, RF, lambda, k)

    %Data från uppgiften
    deltaT_0 = [0;0];
    waterSpecificHeatCapacity = 4186;
    waterDensity = 1020;
    BoxOneEffectiveDepth = 50;
    BoxTwoEffectiveDepth = 2000;
    seconsInAYear = 60*60*24*365; %31 536 000 sekunder
    
    C_1 = (waterSpecificHeatCapacity*BoxOneEffectiveDepth*waterDensity)/seconsInAYear;
    C_2 = (waterSpecificHeatCapacity*BoxTwoEffectiveDepth*waterDensity)/seconsInAYear;
    
    h = 1; %år
    deltaT = zeros(2,length(years)-1);
    deltaT = [deltaT_0, deltaT];
    for t = 1:(length(years)-1)
        deltaT(1,t+h) = deltaT(1,t) + ((RF(t) - (deltaT(1,t))/(lambda) - k*(deltaT(1,t) - deltaT(2,t)))/C_1)*h;
        deltaT(2,t+h) = deltaT(2,t) + ((k*(deltaT(1,t) - deltaT(2,t)))/C_2)*h;
    end
    
end


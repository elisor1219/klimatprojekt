function deltaT = vattenModell(years, RF, lambda, k, C_1, C_2)

    %Data från uppgiften
    deltaT_0 = [0;0];
    
    h = 1; %år
    deltaT = zeros(2,length(years)-1);
    deltaT = [deltaT_0, deltaT];
    for t = 1:(length(years)-1)
        deltaT(1,t+h) = deltaT(1,t) + ((RF(t) - (deltaT(1,t))/(lambda) - k*(deltaT(1,t) - deltaT(2,t)))/C_1)*h;
        deltaT(2,t+h) = deltaT(2,t) + ((k*(deltaT(1,t) - deltaT(2,t)))/C_2)*h;
    end
    
end


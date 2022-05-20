function [B, temperaturCase] = case1(Year, beta, lamda, k, s)
    
    run('utslappRCP45.m')

    U = CO2Emissions;
    %Ändrar år 2022 till 2070
    U(258:306) = linspace(U(258),0,49);
    %Sätter år 2070 - 2200 till noll
    U(306:end) = 0;
    
    [B, temperaturCase] = modellTemp(Year, beta, lamda, k, s, U);

end


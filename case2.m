function [B, temperaturCase] = case2(Year, beta, lamda, k, s)
    
    run('utslappRCP45.m')

    U = CO2Emissions;
    %Ändrar år 2023 och frammåt
    U(259:end) = U(258);
    
    [B, temperaturCase] = modellTemp(Year, beta, lamda, k, s, U);

end
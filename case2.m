function [B, temperaturCase] = case2(Year, beta, lamda, k, s)
    
    run('utslappRCP45.m')

    U = CO2Emissions;
    %�ndrar �r 2023 och framm�t
    U(259:end) = U(258);
    
    [B, temperaturCase] = modellTemp(Year, beta, lamda, k, s, U);

end
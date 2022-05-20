function [B, temperaturCase, C_1, C_2] = case3(Year, beta, lamda, k, s, RF_geoEngineering)
    
    run('utslappRCP45.m')

    U = CO2Emissions;
    %Ändrar år 2023 och fram till och med 2100
    U(259:336) = linspace(U(259),U(259)*2.5,336-258);
    %Ändrar 2101 och frammåt
    U(337:end) = U(259)*2.5;

    if nargin < 6
            [B, temperaturCase, C_1, C_2] = modellTemp(Year, beta, lamda, k, s, U);
    else
            [B, temperaturCase, C_1, C_2] = modellTemp(Year, beta, lamda, k, s, U, RF_geoEngineering);
    end

end
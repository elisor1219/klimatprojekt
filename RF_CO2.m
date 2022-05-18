function outputArg = RF_CO2(year, consCO2)

    run('koncentrationerRCP45.m');
    
    consCO2_noll = CO2ConcRCP45(1)*0.469;

    outputArg = zeros(length(year),1);
    for i = 2:length(year)
        outputArg(i) = 5.35 * log(consCO2(i)/consCO2_noll);
    end

end


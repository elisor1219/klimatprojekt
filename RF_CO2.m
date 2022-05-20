function outputArg = RF_CO2(consCO2)

    run('koncentrationerRCP45.m');
    
    consCO2_noll = CO2ConcRCP45(1);

    outputArg = 5.35 * log(consCO2./consCO2_noll);

end


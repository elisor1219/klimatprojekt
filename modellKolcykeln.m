function B = modellKolcykeln(years, beta, CO2Emissions)
    %st�ller upp alpha, allts� fl�deskoefficienten
    B_noll = [600; 600; 1500]; %GtC
    F_noll = [0, 60, 0
              15, 0, 45
              45, 0, 0]; %GtC/�r
    alpha = F_noll./B_noll;

    %St�ller upp NPP, allts� fl�det mellan box B1 och box B2
    NPP_noll = F_noll(1,2); %Inte helt s�ker p� denna
    NPP = @(B1) NPP_noll * (1 + beta * log(B1/B_noll(1)));

    %St�ller upp utsl�ppet
    U = CO2Emissions;

    %St�ller upp fl�desf�r�ndringen mellan boxaran
    B1_prim = @(t,B1,B2,B3) alpha(3,1)*B3 + alpha(2,1)*B2 - NPP(B1) + U(t);
    B2_prim = @(t,B1,B2,B3) NPP(B1) - alpha(2,3)*B2 - alpha(2,1)*B2;
    B3_prim = @(t,B1,B2,B3) alpha(2,3)*B2 - alpha(3,1)*B3;

    %St�ller up modellen med euler fram�t
    h = 1; %�r
    B = zeros(3,length(years)-1);
    B = [B_noll, B];
    for t=1:(length(years)-1)
        B(1,t+h) =  B(1,t) + B1_prim(t,B(1,t),B(2,t),B(3,t))*h;
        B(2,t+h) =  B(2,t) + B2_prim(t,B(1,t),B(2,t),B(3,t))*h;
        B(3,t+h) =  B(3,t) + B3_prim(t,B(1,t),B(2,t),B(3,t))*h;
    end

end
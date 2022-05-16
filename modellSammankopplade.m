function newModell = modellSammankopplade(years, beta)
    run('utslappRCP45.m')
    
    %B�rjar med att st�lla upp den modell 1 -%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
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
    
    %Fors�tter med att sammanfoga b�da modellerna -%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    %Data fr�n uppgiften
    M_noll = B_noll(1); %GtC
    
    M = zeros(length(years)-1,1);
    M = [M_noll; M];
    h = 1; %�r
    B = zeros(3,length(years)-1);
    B = [B_noll, B];
    antropogenaUtsl = zeros(length(years),1);
    %Hela loppen som modellen simulerar i
    for t = 1:length(years)-1
        antropogenaUtsl(t) = B1_prim(t,B(1,t),B(2,t),B(3,t));
        
        temp_impulssvaret = 0;
        for t_index = 1:t
            temp_impulssvaret = temp_impulssvaret + impulssvaret(t-t_index, t, antropogenaUtsl)*antropogenaUtsl(t_index);
        end
        
        M(t+1) = M(1) + temp_impulssvaret;
        
        %Euler fram�t f�r B_2 och B_3.
        B(1,t+h) = M(t+1);
        B(2,t+h) = B(2,t) + B2_prim(t,B(1,t),B(2,t),B(3,t))*h;
        B(3,t+h) = B(3,t) + B3_prim(t,B(1,t),B(2,t),B(3,t))*h;
        B(4,t+h) = sum(U(1:t)) + sum(B(:,1)) - sum(B(:,t));
    end
    
    newModell = B;
end
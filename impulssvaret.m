function I = impulssvaret(t,t_hatt,U)

    %Data från uppgiften
    A = [0.113, 0.213, 0.258, 0.273, 0.143];
    tau_noll = [2.0, 12.2, 50.4, 243.3, inf];
    k = 3.06E-3;    %Gtc
    
    %Löser summan i (7a) först
    %temp_U = sum(U(1:t));
    temp_U = sum(U(1:t_hatt));
    
    %Löser sedan hela (7a)
    tau = tau_noll .* (1 + k*temp_U);
    
    %Löser varje term i summa för sig
    temp_I = A.*exp(-t./tau);
    
    %Summerar alla termer
    I = sum(temp_I);
    
%     n = length(A);
%     %Alternativt sätt att lösa på%%
%     %Formulerar tidskostnaden
%     tau = @(i,t) tau_noll(i) * (1 + k*sum(U));
%     
%     %Formulerar impulssvaret
%     I_fun = @(t) sum(A(1:n).*exp(-t./tau(1:n,t)));
%     I = I_fun(t);
end
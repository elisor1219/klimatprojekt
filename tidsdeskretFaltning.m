function M = tidsdeskretFaltning(years, U)

    %Data från uppgiften
    M_noll = 600; %GtC
    
    M = zeros(length(years)-1,1);
    for t = 1:length(years)-1
        temp_impulssvaret = 0;
        for t_index = 1:t
            temp_impulssvaret = temp_impulssvaret + impulssvaret(t-t_index, t, U)*U(t_index);
        end
        
        M(t) = M_noll + temp_impulssvaret;
    end
    M = [M_noll; M];
end
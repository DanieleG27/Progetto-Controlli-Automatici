function result = momentoDiInerzia(theta, J_0, J, psi)
    % theta : angolo per il calcolo del momento
    % J : vettore dei coefficienti J_i
    % psi : vettore degli offset angolari ψ_i
    
    sum = J_0;
    for i = 1:length(J)
        sum = sum + J(i) * cos(i * theta + psi(i));
    end
    result = sum;
end



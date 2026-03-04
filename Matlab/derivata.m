function result = derivata(theta, J, psi)
    % theta : angolo per il calcolo del momento
    % J : vettore dei coefficienti J_i
    % psi : vettore degli offset angolari ψ_i
    
    sum = 0;
    for i = 1:length(J)
        sum = sum + J(i)*i * sin(i * theta + psi(i));
    end
    result = -sum;
end
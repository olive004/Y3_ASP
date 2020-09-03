% Q1.2 b) 
function v=rp2(M,N)
    Ar=rand(M,1)*ones(1,N);         % Ar = rand_var_a    U[0,1]     
    Mr=rand(M,1)*ones(1,N);         % Mr = rand_var_b    U[0,1]     
    v=(rand(M,N)-0.5).*Mr+Ar;       % v = rand_var_b(rand(n)-0.5) + rand_var_a 

    % v_2 = rand_var_b w[n] + rand_var_a
    
    % Theo mean & std
    
    
    
    
    
    
    
    
    % WRONK: ASSUMED Ar and Mr WERE rand(n)
    % mean_v = E(v) = E(rand(n)^2 - rand(n)*0.5 + rand(n)) 
    % mean_v = E(rand(n)^2) - 0.5*mean_rand + mean_rand
    % mean_v = integ(x^2 dx)|0,1 - 0.5*0.5 + 0.5
    % mean_v = 1/3 + 0.25 = 1/3+1/4 = 7/12
    % 
    % var = E((v-E(v))^2) = E(v^2) - [E(v)]^2
    % var = E( (rand(n)^2 + 0.5rand(n))^2 ) - (7/12)^2
    % var = E( rand(n)^4 + rand(n)^3 + 0.25*rand(n)^2 ) - (7/12)^2
    % var = E( rand(n)^4 ) + E(rand(n)^3) + 0.25 E(rand(n)^2) - (7/12)^2
    % var = integ(x^4)|0,1 + integ(x^3)|0,1 + 0.25 integ(x^2)|0,1 - (7/12)^2
    % var = 1/5 + 1/4 + 0.25/3 - 49/144
    % var = 0.19305


end
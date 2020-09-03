% Q1.2 c) 
function v=rp3(M,N)
a=0.5;                          
m=3;
v=(rand(M,N)-0.5)*m + a;        % v = 3 (rand(n) - 0.5) + 0.5 = 3rand(n) + 0.5

% v_3 = m w[n] + a

% Theo mean & std
% mean_v = E(v) = E(3 rand(n)) -1.5 + 0.5
% mean_v = 3 mean_rand - 1
% mean_v = 0.5
% 
% var = E((v-E(v))^2) = E(v^2) - [E(v)]^2
% var = E( 9(rand(n))^2 - 6rand(n) + 1) - 0.5^2
% var = 9 integ(x^2 dx)|0,1  - 6 mean_v + 0.25
% var = 9 (x^3/3)|0,1  - 6 * 0.5 + 0.25
% var = 9 (1/3)  - 2.75 = 3 - 2.5 = 0.5

end







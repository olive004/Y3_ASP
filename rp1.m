% Q1.2 a) 
function v=rp1(M,N)
a=0.02;
b=5;
Mc=ones(M,1)*b*sin((1:N)*pi/N);     % Mc = 5 sin( pi/ N)
Ac=a*ones(M,1)*[1:N];               % Ac = 0.02 * n
v=(rand(M,N)-0.5).*Mc+Ac;           % v = 5 sin( n pi/ N) (rand(n) - 0.5) + 0.02 * n

% v = b sin(n pi/ N) (w[n] - 0.5) + a n 


% Theo mean & std
% m_v = E(v) = E(5 sin( pi/ N) (rand(n) - 0.5) + 0.02) = 5*sin(pi/N) E(rand(n)-0.5) +0.02
% m_v = 5*sin(pi/N) E(rand(n)) - 2.5*sin(pi/N) +0.02
% m_v = 5*sin(pi/N) mean_rand - 2.5*sin(pi/N) +0.02
% m_v = 5*sin(pi/N) 0.5 - 2.5*sin(pi/N) +0.02
% m_v = a * n
% 
% FUCK STD ON MATLAB GO TO PAPER
% std_v = E((v-E(v))^2) = E(v^2) - [E(v)]^2
% std_v = E(( b sin( pi/ N) (rand(n) - 0.5) + a )^2) - m_v^2

% std_v = E( b^2(sin( n pi/ N))^2 * (rand(n) - 0.5)^2 + 2a b sin( n pi/ N)
% (rand(n) - 0.5) + a^2) - (a N)^2

% std_v = b^2(sin( n pi/ N))^2 * E((rand(n) - 0.5)^2 + 2a b sin( n pi/ N) (rand(n) - 0.5))
% std_v = 25(sin( pi/ N))^2 * E( rand(n)^2 - rand(n) + 0.25 + 0.2 sin( pi/ N) (rand(n))) - 0.01 sin( pi/ N)
% std_v = 25(sin( pi/ N))^2 * E( rand(n)^2 ) - E(rand(n)) + 0.25 + 0.2 sin( pi/ N)*E(rand(n)) - 0.01 sin( pi/ N)
% std_v = 25(sin( pi/ N))^2 * 1/3 - 0.5 + 0.25 + 0.2 sin( pi/ N)*0.5 - 0.01 sin( pi/ N)
% std_v = 25/3 * (sin( pi/ N))^2 + 0.01 sin( pi/ N) - 0.01 sin( pi/ N) - 0.25
% std_v = 25/3 * (sin( pi/ N))^2 - 0.25

% Joe and Lapo's deriv: 1/12 b^2 (sin(pi n / N))^2 



end
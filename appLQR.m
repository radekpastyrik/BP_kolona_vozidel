function [A,B,C,D,Az,K,H,H_zpet] = appLQR(N)
% clc; clear all; close all;
% N = 3;
pocet_stavu = 2*N-1;
% conds = [0 , 4 , 0 , 2 , 0 , 1 , 0 , 0.4 , 0];
% conds = [0 , 4 , 0 , 2 , 0];
% conds0 = zeros(1,pocet_stavu);

% automaticky matice A,B
A = zeros(2*N-1);
B = zeros(2*N-1,N);

for i = 1:2:(2*N-1)
   A(i,i) = -1;
   A(i+1,i) = 1;
   A(i+1,i+2) = -1;
end
for i = 1:2:(2*N-1)
    B(i,ceil(i/2)) = 1;
end
A = A(1:(2*N-1),1:(2*N-1));
C = eye(pocet_stavu);
% Cy = C;
% Cw = C;
% Cy(2:2:pocet_stavu-1,2:2:pocet_stavu-1) = 0;
% Cw(1:2:pocet_stavu,1:2:pocet_stavu) = 0;
s = tf('s');
D = zeros(pocet_stavu,N);
prenosy =(C*(s*eye(pocet_stavu)-A)^(-1)*B);
Q = eye(pocet_stavu);
Q(2,2) = 4;
Q(3,3) = 2;
Q(4,4) = 4;
Q(5,5) = 2;
if N > 3
    Q(6,6) = 4;
    Q(7,7) = 2;
    if N == 5
        Q(8,8) = 4;
        Q(9,9) = 2;
    end
end

R = eye(N);
[K,~,~] = lqr(A,B,Q,R);
% KKK5,

Az = A-B*K;
% eig(ABK5)
porucha = B; 

% prenosy pro aktualni Q a R
% s = tf('s');
Fc = simplify(C*((s*eye(pocet_stavu)-Az)^(-1))*porucha + D);

H = zeros(1,N-2);    
a=1:2*N-1;
sude=a(mod(a,2)==0);    % sude prenosy odpovidaji z poruchy na odchylky poloh
% podminka retezove stability
for i=1:N-2       % vhodny rozsah   
H(i) = norm(simplify((Fc(sude(i+1),i)/Fc(sude(i),i))),Inf);  % urceni jednotlivych podminek
end
for i=3:N       % vhodny rozsah   
H_zpet(i-2) = norm((Fc(sude(i-2),i)/Fc(sude(i-1),i)),Inf);  % urceni jednotlivych podminek
end
end
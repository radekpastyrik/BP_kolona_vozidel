function [A,B,C,D,Az,K,H,H_zpet] = appSYM(N)
% clear all;close all;clc;
% N = 3;      % pocet vozidel
pocet_stavu = 2*N-1;    % vypocte poce stavu v modelu
% conds = [0 , 4 , 0 , 2 , 0];    % poc. podminky pro 3 vozidla
% if N == 5
%     conds = [0 , 4 , 0 , 3 , 0, 2, 0, -1, 0];
% end

% conds0 = zeros(1,2*N-1);    % nulove poc. podminky

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
Cy = C;
Cw = C;
Cy(2:2:pocet_stavu-1,2:2:pocet_stavu-1) = 0;
Cw(1:2:pocet_stavu,1:2:pocet_stavu) = 0;
D = zeros(pocet_stavu,N);

% pro symetrickou zpetnou vazbu zarucene i prenos z chyby na odchylku rychlost max <= 1 a splnene odchylky polohy
% symetricke rizeni z maple 3 aut
if N == 3
    % dalsi mozna navrzena symetricka reseni pro 3 vozidla
% K = [-8/3 -8/3 1/3 -4/3 1/3;
%     1/3 4/3 -8/3 -4/3 1/3;
%     1/3 4/3 1/3 8/3 -8/3];
% K = [-1 -2/3 0 -1/3 0;
%     0 1/3 -1 -1/3 0;
%     0 1/3 0 2/3 -1];
K = [[0,-1/6,0,-1/12,0];[0,1/12,0,-1/12,0];[0,1/12,0,1/6,0]];
end

% symetricke rizeni z maple 4 aut
if N == 4
    % dalsi mozna navrzena symetricka reseni pro 4 vozidla
% 3 lambdy, 1x ni - temer nulova retezova stab - symetricka zpetna vazba zarucuje temer nulovy hodnoty retez. stability
% K = [[-9/4,-3,3/4,-2,3/4,-1,3/4];[3/4,1,-9/4,-2,3/4,-1,3/4];[3/4,1,3/4,2,-9/4,-1,3/4];[3/4,1,3/4,2,3/4,3,-9/4]];
% K = [[-1,-3/4,0,-1/2,0,-1/4,0];[0,1/4,-1,-1/2,0,-1/4,0];[0,1/4,0,1/2,-1,-1/4,0];[0,1/4,0,1/2,0,3/4,-1]];
K = [[0,-3/16,0,-1/8,0,-1/16,0];[0,1/16,0,-1/8,0,-1/16,0];[0,1/16,0,1/8,0,-1/16,0];[0,1/16,0,1/8,0,3/16,0]];
end

% symetricke rizeni z maple 5 aut
if N == 5
    % dalsi mozna navrzena symetricka reseni pro 5 vozidla
% K = [[-1,-4/5,0,-3/5,0,-2/5,0,-1/5,0];
%     [0,1/5,-1,-3/5,0,-2/5,0,-1/5,0];
%     [0,1/5,0,2/5,-1,-2/5,0,-1/5,0];
%     [0,1/5,0,2/5,0,3/5,-1,-1/5,0];
%     [0,1/5,0,2/5,0,3/5,0,4/5,-1]];

% K = [[-7/10,-4/5,3/10,-3/5,3/10,-2/5,3/10,-1/5,3/10];
%     [3/10,1/5,-7/10,-3/5,3/10,-2/5,3/10,-1/5,3/10];
%     [3/10,1/5,3/10,2/5,-7/10,-2/5,3/10,-1/5,3/10];
%     [3/10,1/5,3/10,2/5,3/10,3/5,-7/10,-1/5,3/10];
%     [3/10,1/5,3/10,2/5,3/10,3/5,3/10,4/5,-7/10]];

K = [[0,-1/5,0,-3/20,0,-1/10,0,-1/20,0];
    [0,1/20,0,-3/20,0,-1/10,0,-1/20,0];
    [0,1/20,0,1/10,0,-1/10,0,-1/20,0];
    [0,1/20,0,1/10,0,3/20,0,-1/20,0];
    [0,1/20,0,1/10,0,3/20,0,1/5,0]];
end

Az = A+B*K;     % dynamika rizeneho systemu
eig(Az);         % jeji vlastni cisla
porucha = B; 

s = tf('s');
Fc = simplify(C*((s*eye(pocet_stavu)-Az)^(-1))*porucha + D);    % prenosy rizeneho sys s poruchou
 
% podminka retezove stability
a=1:2*N-1;
sude=a(mod(a,2)==0);    % sude prenosy odpovidaji z poruchy na odchylky poloh
for i=1:N-2       % vhodny rozsah   
H(i) = norm((Fc(sude(i+1),i)/Fc(sude(i),i)),Inf);  % urceni jednotlivych podminek
end

for i=3:N       % vhodny rozsah   
H_zpet(i-2) = norm((Fc(sude(i-2),i)/Fc(sude(i-1),i)),Inf);  % urceni jednotlivych podminek
end
end


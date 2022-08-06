function [A,B,C,D,Az,K,H,H_zpet] = appSZV(N,topologie)
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

% rizeni dle topologie BD,PF,opak PF pro 3 auta
if N == 3
    if topologie == 1
        K = [[-1/2,-1/2,0,0,0];[3/4,3/4,-3/2,0,0];[0,0,1,1,-1]]; % PF
    elseif topologie == 2
        K = [[0,-1/4,7/18,1/18,0];[0,0,-1,-1,7/4];[0,0,0,0,1/2]]; %SF druhy smer-opak
    else 
        K = [[-7/3,-8/3,4/3,-5/2,0];[2/3,4/3,-1/6,1/2,-3/2];[0,0,3/2,1/2,-5/2]]; % BD
    end    
end

if N == 4
    if topologie == 1
%         K = 0; % PF
    elseif topologie == 2
        K = [[-2,-2,1,0,0,0,0];[0,0,-2,-2,1,0,0];[0,0,0,0,-4,-6,3];[0,0,0,0,0,0,-1]]; % SF
    else
        K = [[-5/2,-3,7/4,1/2,0,0,0];[-1/7,-2/7,-7/2,-37/7,18/7,0,0];[0,0,-7/4,-7/2,1,1,-1/2];[0,0,0,0,2,4,-3]]; %BD
    end   
end

if N == 5
    if topologie == 1
%         K = 0; % PF
    elseif topologie == 2
%         K = 0; % SF
    else
%         K = 0; %BD
    end   
end

Az = A+B*K;     % dynamika rizeneho systemu
eig(Az);         % jeji vlastni cisla
porucha = B; 

s = tf('s');
Fc = simplify(C*((s*eye(pocet_stavu)-Az)^(-1))*porucha + D);    % prenosy rizeneho sys s poruchou
 
% podminka retezove stability dopredu
a=1:2*N-1;
sude=a(mod(a,2)==0);    % sude prenosy odpovidaji z poruchy na odchylky poloh
for i=1:N-2       % vhodny rozsah   
H(i) = norm(simplify((Fc(sude(i+1),i)/Fc(sude(i),i))),Inf);  % urceni jednotlivych podminek
end

% podminka retezove stability dozadu - opak dopredne ss
for i=3:N       % vhodny rozsah   
H_zpet(i-2) = norm((Fc(sude(i-2),i)/Fc(sude(i-1),i)),Inf);  % urceni jednotlivych podminek
end

H      %vypis vycislene podminek retezove stability
H_zpet
end
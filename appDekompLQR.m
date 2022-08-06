function [A,B,C,D,Az,K,H,H_zpet] = appDekompLQR(N)
% N = 3;
pocet_stavu = 2*N - 1;
%strukturalni dekompozice - lqr jako distribuovane rizeni
% automaticky matice A,B
A = zeros(pocet_stavu);
B = zeros(pocet_stavu,N);

% conds = [0 , 4 , 0 , 2 , 0];    % poc. podminky pro 3 vozidla
% if N == 5
%     conds = [0 , 4 , 0 , 3 , 0, 2, 0, -1, 0];
% end
% 
% conds0 = zeros(1,pocet_stavu);    % nulove poc. podminky

for i = 1:2:(pocet_stavu)
   A(i,i) = -1;
   A(i+1,i) = 1;
   A(i+1,i+2) = -1;
end
for i = 1:2:(pocet_stavu)
    B(i,ceil(i/2)) = 1;
end
A = A(1:(pocet_stavu),1:(pocet_stavu));

C = eye(pocet_stavu);
Cy = C;
Cw = C;
Cy(2:2:pocet_stavu-1,2:2:pocet_stavu-1) = 0;
Cw(1:2:pocet_stavu,1:2:pocet_stavu) = 0;
D = zeros(pocet_stavu,N);

%leader
if N == 6
    Q_lead = 1;
    R_lead = 1; 
elseif N == 7
    Q_lead = 5;
    R_lead = 2; 
else
    Q_lead = 10;
    R_lead = 2; 
end

a_lead = A(1,1);
b_lead = B(1,1);

[K_lead,~,~] = lqr(a_lead,b_lead,Q_lead,R_lead);
K_lead1 = [K_lead,0,0];
sys0 = ss(a_lead,b_lead,1,0);
Az_lead = a_lead - b_lead*K_lead1;

K = zeros(N,pocet_stavu);
K(1,1) = K_lead;
ai = A(1:3,1:3);
bi = B(1:3,1);
Bi = B(3:5,3);
Qi = eye(3);
if N == 5
    Qi(2,2) = 20;
    Ri = 2; 
elseif N == 6
	Qi(2,2) = 1;
	Ri = 1; 
elseif N == 7
	Qi(2,2) = 25;
	Ri = 3; 
else
    Qi(2,2) = 20;
    Ri = 3;  
end

odkud = 1;
kam = odkud+2;
for i=2:N
    if i == 2
      Abk = ai - bi*K(i-1,odkud:kam) ;   
    else
      Abk = ai - bi*K(i-1,odkud-2:kam-2); 
    end
    [K1o,~,~] = lqr(Abk,Bi,Qi,Ri);
    K(i,odkud:kam) = K1o(:);
    odkud = kam;
    kam = kam+2;
end

Az = A-B*K;     % dynamika rizeneho systemu
ei = eig(Az);        % jeji vlastni cisla
porucha = B; 

s = tf('s');
Fc = simplify(C*((s*eye(pocet_stavu)-Az)^(-1))*porucha + D);    % prenosy rizeneho sys s poruchou
 
% podminka retezove stability
a=1:2*N-1;
sude=a(mod(a,2)==0);    % sude prenosy odpovidaji z poruchy na odchylky poloh
for i=1:N-2       % vhodny rozsah   
H(i) = norm((Fc(sude(i+1),i)/Fc(sude(i),i)),Inf);  % urceni jednotlivych podminek
end
% podminka retezove stability dozadu -vhodnejsi- opak dopredne ss
for i=3:N       % vhodny rozsah   
H_zpet(i-2) = norm((Fc(sude(i-2),i)/Fc(sude(i-1),i)),Inf);  % urceni jednotlivych podminek
end
end

                 

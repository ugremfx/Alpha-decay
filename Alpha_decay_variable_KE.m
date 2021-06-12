clear all
clc
close all

no_seg = 10; % Number of segments the barrier is divided in Change this

R0 = 1.4; % fm , given in assignment comments Change this
AD = 234; % Mass number of Daughter nucleus Change this
R_in = R0*AD^(1/3); % Radius of the potential well Final radius of the potential well

mass_alpha = 3727.3794066; % MeV/c^2 (megaelectronvolts)
mass_daughter = 218010.234 ; % MeV/c^2 Change this
mass_parent = 238.05078826*931.4941024228; % MeV/c^2



Z1 = 92; % Daughter nucleus. Change this
Z2 = 2;  % alpha particle

% Daughter nucleus: Thorium
fine_str_ct = 1/137.035399; % [-]
hbarc = 197.3269631; % MeV*fm
% hbar = 6.5821195695091*10^-22; % [MeV s] (Megaelectronvolt seconds) Plank's constant

% KE_alpha = 4.19; % MeV;  Source: http://physics-database.group.shef.ac.uk/phy303/phy303-4.html
% KE_alpha = (-1801.68856 + 28.296 + 1777.662666)*234/238; % other notes
E_kin = mass_parent-mass_daughter-mass_alpha; % Decay energy from mass balance
KE_alpha = E_kin; % Assumption: The kinetic energy of alpha is the same as the energy decay. The α-particle emerges with a kinetic energy T_alpha which is slightly below the value of Q.  This is because if the parent nucleus is at rest before decay there must be some recoil of the daughter nucleus in order to conserve momentum.

M_mu = (mass_alpha*mass_daughter)/(mass_daughter+mass_alpha); % [MeV/c^2]
v_rel = sqrt(2*E_kin/M_mu); % [c]
v_rel = v_rel*2.99792458e23;  % Transform to femtometers/s

R_out = fine_str_ct*Z1*Z2*hbarc/KE_alpha; % This will be in fm. It represents the outer limit of quantum tunneling (outer limit of barrier penetration)

x_barrier = R_in:((R_out-R_in)/(no_seg)):R_out; % x positions of the barrier interfaces
len_seg = x_barrier(2)-x_barrier(1);

V0 = 134; % MeV % See the drawing in alphanotes.pdf
Veff = [-V0+fine_str_ct*Z1*Z2*hbarc/R_in , fine_str_ct*Z1*Z2*hbarc./(x_barrier(1:end-1) + len_seg/2), 0]; %    x_barrier

%figure
%plot([x_barrier,80],Veff,'*')

k(1) =  sqrt(2*mass_alpha*(E_kin-2*Veff(1)))/hbarc; % [fm^-1] This is the first k 
k(2:no_seg+1) = sqrt(2*mass_alpha*(2*Veff(2:end-1)-E_kin))/hbarc; % [fm^-1] These are the k's in between
k(no_seg+2) =  sqrt(2*mass_alpha*(E_kin-2*Veff(end)))/hbarc; % [fm^-1] This is the last k
%Note: Number of k's = Number of segments + 2

A = zeros((2*no_seg+4),(2*no_seg+4)); % Pre allocation for AX=B This is a matrix populated with zeros

for ii = 2:no_seg %Writes all the equations for segments except for the first and last segment
    %These eight lines populate the matrix with the coefficients
    A(2*ii, 2*ii-1) = exp(k(ii)*x_barrier(ii)); %Writes the SE for iith barrier
    A(2*ii, 2*ii  ) = exp(-k(ii)*x_barrier(ii)); 
    A(2*ii, 2*ii+1) = -exp(k(ii+1)*x_barrier(ii));
    A(2*ii, 2*ii+2) = -exp(-k(ii+1)*x_barrier(ii));
    
    A(2*ii+1, 2*ii-1) = k(ii)*exp(k(ii)*x_barrier(ii)); %Writes the derivative of SE for iith barrier
    A(2*ii+1, 2*ii  ) = -k(ii)*exp(-k(ii)*x_barrier(ii));
    A(2*ii+1, 2*ii+1) = -k(ii+1)*exp(k(ii+1)*x_barrier(ii));
    A(2*ii+1, 2*ii+2) =  k(ii+1)*exp(-k(ii+1)*x_barrier(ii));
end
A(1,1) = 1; %Sets A = 1
A(end,end) = 1; %This is G. The incoming wave from the rightmost barrier, which is equal to zero

A(2, 1) = exp(1i*k(1)*x_barrier(1)); %These four lines populate the matrix
A(2, 2) = exp(-1i*k(1)*x_barrier(1));%with the SE for the first segment of the 
A(2, 3) = -exp(k(2)*x_barrier(1));%barrier.
A(2, 4) = -exp(-k(2)*x_barrier(1));%Eqn. 13, 14 and 15

A(3, 1) = 1i*k(1)*exp(1i*k(1)*x_barrier(1));%These four lines populate the matrix 
A(3, 2) = -1i*k(1)*exp(-1i*k(1)*x_barrier(1));%with the derivative of the SE
A(3, 3) = -k(2)*exp(k(2)*x_barrier(1));%for the first segment of the 
A(3, 4) =  k(2)*exp(-k(2)*x_barrier(1));%barrier. %Eqn. 13, 14 and 15

ii = (no_seg+1);
A(2*(no_seg+1), 2*(no_seg+1)-1) = exp(k(ii)*x_barrier(ii));%These lines populate the matrix
A(2*(no_seg+1), 2*(no_seg+1)  ) = exp(-k(ii)*x_barrier(ii));%with the SE
A(2*(no_seg+1), 2*(no_seg+1)+1) = -exp(1i*k(ii+1)*x_barrier(ii));%for the 
A(2*(no_seg+1), 2*(no_seg+1)+2) = -exp(-1i*k(ii+1)*x_barrier(ii));%last segment

A(2*(no_seg+1)+1, 2*(no_seg+1)-1) = k(ii)*exp(k(ii)*x_barrier(ii));%These lines populate the matrix
A(2*(no_seg+1)+1, 2*(no_seg+1)  ) = -k(ii)*exp(-k(ii)*x_barrier(ii));%with the derivative of the SE
A(2*(no_seg+1)+1, 2*(no_seg+1)+1) = -1i*k(ii+1)*exp(1i*k(ii+1)*x_barrier(ii));%for the
A(2*(no_seg+1)+1, 2*(no_seg+1)+2) = 1i*k(ii+1)*exp(-1i*k(ii+1)*x_barrier(ii));%last segment Eqn 13, 14 and 15
    
B = [1;zeros(2*no_seg+3,1)]; %These is the constants of the systems of equation

% Solve

X = A\B; % A B C1 D1 C2 D2 ... F G This solves the system of equations AX = B 

T = (abs(X(end-1)))^2*k(end)/k(1); %Here we calculate the transmission coefficient

tau = (T*v_rel/(2*R_in))^(-1); % [s] Here we calculate tau
half_time = tau*log(2); % [s] Here we calculate the half-time
fprintf('Half time is %6.2e seconds \n',half_time) %Here we print the half-time

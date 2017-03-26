% This MATLAB program checks the feasibility of LMIs from Theorems 1 and 2 of the paper 
% A. Selivanov and E. Fridman, "Predictor-based networked control in the presence of uncertain time-varying delays," in 55th IEEE Conference on Decision and Control, 2016, pp. 501–506.

%% LMIs of Theorem 1
% System parameters
M=10;   % the cart mass
m=1;    % the pendulum mass
l=3;    % the length of the pendulum arm
g=10;   % the gravitational acceleration

A=[0 1 0 0; 0 0 -m*g/M 0; 0 0 0 1; 0 0 g/l 0]; 
B=[0; 1/M; 0; -1/(M*l)]; 
K=[2 12 378 210]; 

h=.0368; r0=.2; etaM=.01; r1=.2; muM=.01; alpha=.01; 
if LMI_CDC16a_th1(A,B,K,h,r0,etaM,r1,muM,alpha)
    display('Feasible'); 
else
    display('Not Feasible'); 
end

%% LMIs of Theorem 2
% System parameters
A=[0 1; 0 -.1]; 
B=[0; .1]; 
C=[1 0]; 
K=-[3.75 11.5]; 
L=-[1.4 .36]';

h=.06; r0=1; etaM=0.1; r1=1; muM=0.1; alpha=.001; sigma=.07; 
display(['Omega=' num2str(LMI_CDC16a_th2(A,B,C,K,L,h,r0,etaM,r1,muM,alpha,sigma))]); 

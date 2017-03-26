function flag=LMI_CDC16a_th1(A,B,K,h,r0,etaM,r1,muM,alpha)
% This MATLAB program checks the feasibility of LMIs from Theorem 1 of the paper 
% A. Selivanov and E. Fridman, "Predictor-based networked control in the presence of uncertain time-varying delays," in 55th IEEE Conference on Decision and Control, 2016, pp. 501–506.

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)
% and SeDuMi solver (http://sedumi.ie.lehigh.edu/)

% Input: 
% A,B   - the parameters of the system (5); 
% K     - the controller gain from (9); 
% h     - the sampling period; 
% r0    - the sensors-to-controller known constant delay; 
% etaM  - a bound for the sensors-to-controller time-varying delay uncertainty; 
% r1    - the controller-to-actuators known constant delay; 
% muM   - a bound for the controller-to-actuators time-varying delay uncertainty; 
% alpha - a desired decay rate; 

% Output: 
% flag==1 if LMIs are feasible, flag==0 otherwise. 

n=size(A,1); 

%% Decision variables 
P=sdpvar(n); 
S=sdpvar(n); 
S0=sdpvar(n); 
S1=sdpvar(n); 
R0=sdpvar(n); 
R1=sdpvar(n); 
P2=sdpvar(n,n,'f'); 
P3=sdpvar(n,n,'f'); 
G0=sdpvar(n,n,'f'); 
G1=sdpvar(n,n,'f'); 
G2=sdpvar(n,n,'f'); 
G3=sdpvar(n,n,'f'); 

%% Notations
tauBar=h+etaM; 
tauM=r0+r1+h+etaM+muM; 
rhoBar=exp(-2*alpha*tauBar); 
rhoM=exp(-2*alpha*tauM); 

%% The main LMI (Phi<=0)
Phi=blkvar; 
Phi(1,1)=2*alpha*P+S0-rhoBar*R0+P2'*A+A'*P2; 
Phi(1,2)=P-P2'+A'*P3; 
Phi(1,3)=rhoBar*(R0-G0)+P2'*B*K; 
Phi(1,4)=rhoBar*G0; 
Phi(1,6)=-P2'*expm(A*(r0+r1))*B*K; 
Phi(1,7)=P2'*expm(A*(r0+r1))*B*K; 
Phi(2,2)=tauBar^2*R0+(tauM-r0-r1)^2*R1-P3-P3'; 
Phi(2,3)=P3'*B*K; 
Phi(2,6)=-P3'*expm(A*(r0+r1))*B*K; 
Phi(2,7)=P3'*expm(A*(r0+r1))*B*K; 
Phi(3,3)=-rhoBar*(2*R0-G0-G0'); 
Phi(3,4)=rhoBar*(R0-G0); 
Phi(4,4)=rhoBar*(S-S0-R0); 
Phi(5,5)=exp(-2*alpha*(r0+r1))*(S1-S)-rhoM*R1; 
Phi(5,6)=rhoM*(R1-G1); 
Phi(5,7)=rhoM*(G1-G2); 
Phi(5,8)=rhoM*G2; 
Phi(6,6)=-rhoM*(2*R1-G1-G1'); 
Phi(6,7)=rhoM*(R1-G1+G2-G3); 
Phi(6,8)=rhoM*(G3-G2); 
Phi(7,7)=-rhoM*(2*R1-G3-G3'); 
Phi(7,8)=rhoM*(R1-G3); 
Phi(8,8)=-rhoM*(S1+R1); 
Phi=sdpvar(Phi); 

%% Park's conditions 
Park0=[R0 G0; G0' R0]; 
Park1=[R1 G1; G1' R1]; 
Park2=[R1 G2; G2' R1]; 
Park3=[R1 G3; G3' R1]; 

%% Solution of LMIs
LMIs=[Phi<=0, P>=0, S>=0, S0>=0, S1>=0, R0>=0, R1>=0, Park0>=0, Park1>=0, Park2>=0, Park3>=0]; 
options=sdpsettings('solver','sedumi','verbose',0);
sol=optimize(LMIs,[],options); 

flag=0; 
if sol.problem==0
    [primal,~]=check(LMIs); % Checking that the solver returned a proper solution
    flag=(min(primal)>=0 && primal(2)>0); 
else
    yalmiperror(sol.problem) 
end

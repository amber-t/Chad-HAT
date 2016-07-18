%% paramters (days)
%humans
H=14500; %human population size
muH=4.66e-5; %human death rate
betaH=0.1751; %transmission prob from tsetse > human
taoH=1/12; %human incubation period
gammaH1=1/526; %stage 1 infectious period
gammaH2=1/252; %stage 2 infectious period

%tsetse
V=17*14500; %tsetse population size
eta=1/20; %pupae stage duration
Bv=0.05; %tsetse birth rate
muV0=0.030; %tsetse death rate no comp
muV1=0.0002; %tsetse death comp param
sigmaV=1; %susceptibility period
a=0.333; %human bite rate
betaVH=0.3750; %proportion of tsetse bites on humans
taoV=1/25; %incubation period in tsetse
betaV=0.2; %transmission prob from humans to tsetse

%treatment
epi=0.965; %efficiency of stage II treatment (nifurtimox-eflornithine)
zeta=1/202; %treatment seeking rate for stage II patietns
p=0.007; %probability of death due to stage 2 treatment failure
delta=1/50; %immune period in humans after treatment
rho=0.87; %sensitivity of CATT-TL


%% initial conditions
%humans
X_0=zeros(10,1);
X_0(1)=0.75; %suscep
X_0(2)=0.2; %exposed
X_0(3)=0.05; %infectious s1
X_0(4)=0; %infectious s2
X_0(5)=0; %recovered

X_0 = H*X_0;

%tsetse
X_0(6)=0.1; %pupae stage
X_0(7)=0.5; %susceptible
X_0(8)=0.1; %exposed
X_0(9)=0; %infected
X_0(10)=0.3; %recovered




%% ode45
[t,X]=ode15s(@HAT,[0:1:500],X_0,[],delta,a,betaVH,betaH,H,muH,taoH,gammaH1,epi,zeta,...
    gammaH2,p,Bv,eta,sigmaV,taoV,rho,muV0,muV1,V,betaV);

Hs=X(:,1);
He=X(:,2);
Hi1=X(:,3);
Hi2=X(:,4);
Hr=X(:,5);

Vp=X(:,6);
Vs=X(:,7);
Ve=X(:,8);
Vi=X(:,9);
Vr=X(:,10);

plot(t,Hs,t,He,t,Hi1,t,Hi2,t,Hr)
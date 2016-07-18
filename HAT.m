%% HAT: SEIR model with vectored transmission
function dX=HAT(t,X,delta,a,betaVH,betaH,H,muH,taoH,gammaH1,epi,zeta,...
    gammaH2,p,Bv,eta,sigmaV,taoV,rho,muV0,muV1,V,betaV)
%% parameters
Hs=X(1);
He=X(2);
Hi1=X(3);
Hi2=X(4);
Hr=X(5);

Vp=X(6);
Vs=X(7);
Ve=X(8);
Vi=X(9);
Vr=X(10);


Bh = muH*H + (1-rho)*gammaH2*Hi2 + rho*(1-epi)*p*zeta*Hi2; %human birth rate
lambdaVH = betaV*betaVH*(Hi1/H); %prob of a suscpetible tsetse become infected from a human blood meal
muV = muV0*(1 + muV1*V); %tsetse death rate

dX=zeros(10,1);

%% differential equations
%humans
dHs = Bh + delta*Hr - a*betaVH*betaH*Vi*(Hs/H) - muH*Hs;
dHe = a*betaVH*betaH*Vi*(Hs/H) - taoH*He - muH*He;
dHi1 = taoH*He - gammaH1*Hi1 - muH*Hi1;
dHi2 = gammaH1*Hi1 - rho*epi*zeta*Hi2 - (1 - rho)*gammaH2*Hi2 - rho*(1-epi)*p*zeta*Hi2 - muH*Hi2;
dHr = rho*epi*zeta*Hi2 - delta*Hr - muH*Hr;

%tsetse
dVp = Bv*V + eta*Vp;
dVs= eta*Vp - a*Vs - sigmaV*Vs - muV*Vs;
dVe = a*lambdaVH *Vs - taoV*Ve - muV*Ve;
dVi = taoV*Ve - muV*Vi;
dVr = a*(1-lambdaVH)*Vs + sigmaV*Vs - muV*Vr;

dX(1)=dHs;
dX(2)=dHe;
dX(3)=dHi1;
dX(4)=dHi2;
dX(5)=dHr;

dX(6)=dVp;
dX(7)=dVs;
dX(8)=dVe;
dX(9)=dVi;
dX(10)=dVr;

%% HAT: SEIR model with vectored transmission
function dX=HAT(t,X,delta,a,betaVH,betaH,betaV,muH,taoH,gammaH1,epi1,epi2,zeta1,zeta2,...
    gammaH2,p2,Bv,eta,sigmaV,taoV,muV0,muV1,P1,P1PD,P1TP,P2,P2PD,P2TP)
%% parameters
Hs=X(1); He=X(2); Hi1=X(3); Hi2=X(4); Hr=X(5); 
Vp=X(6); Vs=X(7); Ve=X(8); Vi=X(9); Vr=X(10);

H=Hs+He+Hi1+Hi2+Hr;
V=Vp+Vs+Ve+Vi+Vr;

phi1=P1*P1PD*P1TP; %coverage stage 1 
phi2=P2*P2PD*P2TP; %coversage stage 2

Bh = muH*H + (1-phi2)*gammaH2*Hi2 + phi2*(1-epi2)*p2*zeta2*Hi2; %human birth rate
lambdaVH = betaV*betaVH*(Hi1/H); %prob of a suscpetible tsetse become infected from a human blood meal
muV = muV0*(1 + muV1*V); %tsetse death rate


dX=zeros(10,1);

%% differential equations
%humans
dHs = Bh + delta*Hr - a*betaVH*betaH*Vi*(Hs/H) - muH*Hs;
dHe = a*betaVH*betaH*Vi*(Hs/H) - taoH*He - muH*He;
dHi1 = taoH*He - phi1*epi1*zeta1*Hi1 - (1-phi1)*gammaH1*Hi1 - muH*Hi1;
dHi2 = (1-phi1)*gammaH1*Hi1 + - phi2*epi2*zeta2*Hi2 - (1 - phi2)*gammaH2*Hi2 - phi2*(1-epi2)*p2*zeta2*Hi2 - muH*Hi2;
dHr = phi1*epi1*zeta1*Hi1 + phi2*epi2*zeta2*Hi2 - delta*Hr - muH*Hr;

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

%% HAT: SEIR model with vectored transmission
function dX=HAT(t,X,deltaH,aH,betaVH,betaH,muH,tauH,gammaH1,eps1,eps2,zeta1,zeta2,...
    gammaH2,p2,BV,eta,sigmaV,tauV,muV0,muV1,P1,P1PD,P1TP,P2,P2PD,P2TP,betaV)
%% parameters
Hs=X(1); He=X(2); Hi1=X(3); Hi2=X(4); Hr=X(5); 
Vp=X(6); Vs=X(7); Ve=X(8); Vi=X(9); Vr=X(10);

H=Hs+He+Hi1+Hi2+Hr;
V=Vs+Ve+Vi+Vr;

phi1=P1*P1PD*P1TP; %coverage stage 1 
phi2=P2*P2PD*P2TP; %coversage stage 2

Bh = muH*H + (1-phi2)*gammaH2*Hi2 + phi2*(1-eps2)*p2*zeta2*Hi2; %human birth rate
lambdaVH = betaV*betaVH*(Hi1/H); %prob of a suscpetible tsetse become infected from a human blood meal
muV = muV0*(1 + muV1*V); %tsetse death rate


dX=zeros(10,1);

%% differential equations
%humans
dHs = Bh + deltaH*Hr - aH*betaVH*betaH*Vi*(Hs/H) - muH*Hs;
dHe = aH*betaVH*betaH*Vi*(Hs/H) - tauH*He - muH*He;
dHi1 = tauH*He - phi1*eps1*zeta1*Hi1 - (1-phi1)*gammaH1*Hi1 - muH*Hi1;
dHi2 = (1-phi1)*gammaH1*Hi1 - phi2*eps2*zeta2*Hi2 - (1 - phi2)*gammaH2*Hi2 - phi2*(1-eps2)*p2*zeta2*Hi2 - muH*Hi2;
dHr = phi1*eps1*zeta1*Hi1 + phi2*eps2*zeta2*Hi2 - deltaH*Hr - muH*Hr;

%tsetse
dVp = BV*V - eta*Vp;
dVs= eta*Vp - aH*Vs - sigmaV*Vs - muV*Vs;
dVe = aH*lambdaVH *Vs - tauV*Ve - muV*Ve;
dVi = tauV*Ve - muV*Vi;
dVr = aH*(1-lambdaVH)*Vs + sigmaV*Vs - muV*Vr;

dX = vertcat(dHs,dHe,dHi1,dHi2,dHr,dVp,dVs,dVe,dVi,dVr);

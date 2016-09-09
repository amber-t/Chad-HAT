%% Data and Sample Size
tic;

Data = [273, 139;
       197, 68;
       53 , 18;
       96 , 36]; %[s1, s2]; 2002, 2004, 2006, 2012

SampSize = [11046, 11046;
         13783, 13783;
         9541 , 9541;
         13410, 13410]; %sample size; 2002, 2004, 2006, 2012

x = betarnd(273,11046,10000,1);
y = betarnd(139,11046,10000,1);
ci1 = quantile(x,[0.0025,0.9975]); %s1, 2002 confidence interval
ci2 = quantile(y,[0.0025,0.9975]); %s2, 2002 confidence interval

%% Likelihood

N=1000000;
params=zeros(N,3);

parfor j = 1:N
       betaH = rand;
       betaVH = 0.1+0.5*rand;
       zeta=0.7*rand;

       %run model with randomly chosen parameters
       params(j,:) = [betaH, betaVH, zeta];
       out = HATrun(params(j,:)); %run model

       if (out{1}(1)<=ci1(1)) || (out{1}(1)>=ci1(2)) || (out{1}(5)<=ci2(1))...
               || (out{1}(5)>=ci2(2));
           Likelihood(j)=0;
       else
           Lik1=1; Lik2=1;
           for i=1:4
               Lik1=Lik1*betapdf(out{1}(i),Data(i,1),SampSize(i,1));
               Lik2=Lik2*betapdf(out{1}(i+4),Data(i,2),SampSize(i,2));
               % Lik1(i)=betapdf(out{1}(i),Data(i,1),SampSize(i,1));
               % Lik2(i)=betapdf(out{1}(i+4),Data(i,2),SampSize(i,2))
           end
           Likelihood(j)=Lik1*Lik2;
           %Likelihood(j)=prod(Lik1)*prod(Lik2);
      end
end

save('output4','params','Likelihood')

toc

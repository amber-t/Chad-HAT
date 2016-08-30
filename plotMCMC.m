%% find non-zero likelihood parameters
load('output.mat')
a=find(Likelihood~=0); %index of non-zero likelihoods
b=Likelihood(Likelihood~=0); %non-zero likelihood values
p=params(a,:); %parameter values of non-zero likelihood [betaH,betaVH,zeta]
for i=1:20
out=HATrun(p(i,:));
m=out{1};
s1(:,i)=m(:,1);
s2(:,i)=m(:,2);
end

%% plots
Data = [273, 139;
       197, 68;
       53 , 18;
       96 , 36]; %[s1, s2]; 2002, 2004, 2006, 2012

SampSize = [11046, 11046;
            13783, 13783;
            9541 , 9541;
            13410, 13410]; %sample size; 2002, 2004, 2006, 2012

s1_data = Data(:,1)./SampSize(:,1);
s2_data = Data(:,2)./SampSize(:,2);

%legendCell=cellstr(num2str(b','b=%-d'));

%stage 1
figure
plot([2,4,6,12],s1_data,'o')
hold on
plot([2,4,6,12],s1)

title('Stage 1 Model Fit')
xlabel('year')
ylabel('proportion stage 1')
legend('data','1.1e-06','7.3e-06','7.9e-09','8.0e-01','2.1e+05','1.7e-05','1.4e+04','2.8e+02','1.7e-01','3.1e+03','4.2e+03',...
    '2.8e+01','1.2e+02','2.7e+03','8.7e+01','3.4e+02','2.4e-03','1.0e+03','2.4e+01','2.8e+02','2.3e-06','7.3e+02')

%stage 2
figure
plot([2,4,6,12],s2_data,'o')
hold on
plot([2,4,6,12],s2)

title('Stage 2 Model Fit')
xlabel('year')
ylabel('proportion stage 1')
legend('data','1.1e-06','7.3e-06','7.9e-09','8.0e-01','2.1e+05','1.7e-05','1.4e+04','2.8e+02','1.7e-01','3.1e+03','4.2e+03',...
    '2.8e+01','1.2e+02','2.7e+03','8.7e+01','3.4e+02','2.4e-03','1.0e+03','2.4e+01','2.8e+02','2.3e-06','7.3e+02')

%likelihoods
figure
stem(b)
title('Likelihood')
%% find non-zero likelihood parameters
load('output3.mat')
a=find(Likelihood~=0); %index of non-zero likelihoods
b=Likelihood(Likelihood~=0); %non-zero likelihood values
p=params(a,:); %parameter values of non-zero likelihood [betaH,betaVH,zeta]
s=size(b);
parfor i=1:s(2)
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

x=[2,4,6,12]; %2002,2004,2006,2012
X=[x, flip(x)]; %x vals for fill

%stage 1
subplot(1,2,1)
%plot(x,s1,'color',[0.7 0.7 0.7])
[y_min,i_min]=min(s1(1,:));
[y_max,i_max]=max(s1(1,:));
Y=[s1(:,i_min)', flip(s1(:,i_max)')];
fill(X,Y,[0.75 0.75 0.75],'EdgeColor','none')
hold on

for j=1:4
r1(:,j)=betarnd(Data(j,1),SampSize(j,1),10000,1); %generate 95 CIs
err1(j,:)=quantile(r1(:,j),[0.025,0.975]);
end
errorbar(x,s1_data,err1(:,1),err1(:,2),'ok','linewidth',1.5) %data s1, errorbars

title('Stage 1')
xlabel('year')
ylabel('proportion stage 1')

%stage 2
subplot(1,2,2)
%plot(x,s2,'color',[0.7 0.7 0.7])
[y2_min,i2_min]=min(s2(1,:));
[y2_max,i2_max]=max(s2(1,:));
Y2=[s2(:,i2_min)', flip(s2(:,i2_max)')];
fill(X,Y2,[0.75 0.75 0.75],'EdgeColor','none')
hold on

for j=1:4
r2(:,j)=betarnd(Data(j,2),SampSize(j,2),10000,1); %generate 95 CIs
err2(j,:)=quantile(r2(:,j),[0.025,0.975]);
end
errorbar(x,s2_data,err2(:,1),err2(:,2),'ok','linewidth',1.5) %data s2, errorbars

title('Stage 2')
xlabel('year')
ylabel('proportion stage 1')

%likelihoods
figure
stem(b)
title('Likelihood')


% parfor i = 1:length(B)
%     weights(i) = L(i)/sum(L);
% end

% total = 500;

% j = 1;
% while j < total+1
%     k =randi(length(B),1);
%     if(rand <weights(k))
%     posterior(j,:)=Xpar(k,:);
%      j = j+1;
%     end
% end

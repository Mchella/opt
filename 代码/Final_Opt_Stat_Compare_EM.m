clear
clc
[A,degree] = Get_arenas_email_Network();

if ~exist('Results_OnExcels')
    mkdir('Results_OnExcels')
end

T = 10;
h = 0.01;
xtk = 0:h:T;
n = length(xtk);
maxTimes = 9;
theta_init = 0;

alpha=0.1;
beta=0.1;
c1=750;
c2=1;
l=0;
t=datestr(now,'yyyymmddTHHMMSS');
%最优控制
[ EP,OJ,utheta,OJS] = WAS_OptCont_Euler(T,h,maxTimes,A,theta_init,alpha,beta,c1,c2);
%静态控制
for i=1:1:100
    l=l+1;
    [ staEP,staOJ(l)] = WAS_StaCont(T,h,A,alpha,beta,c1,c2);
end
staOJ=[nan staOJ];

%%
figure;
subplot(1,2,1)
plot(xtk,utheta,'g','linewidth',2);
xlabel('Time');
ylabel('Control');

subplot(1,2,2)
stem(OJ);
hold on;
stem(staOJ);
ylabel('J(\bf{w})');

saveas(gcf,['Fig_OnNet\EM',t,'.jpg']);
% %%
Igor = {101,2};
Igor{1,1}=1;
for i=2:101
    Igor{i,1}=i+5;
end
Igor{1,2} =OJ;
for i=2:101
    Igor{i,2}=staOJ(i);
end
 filename = ['Results_OnExcels\累积效益成本EM',t,'.xls'];
  xlswrite(filename,Igor);
  
  Igor={n-1,2};
  for i=1:n-1
      Igor{i,1} = xtk(i);
      Igor{i,2} = utheta(i);
  end
filename = ['Results_OnExcels\thetaEM',t,'.xls'];
xlswrite(filename,Igor); 
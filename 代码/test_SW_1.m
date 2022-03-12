clear
clc
A = load('SW.txt'); % 小世界网络邻接矩阵导入A
if ~exist('Results_OnExcels')
    mkdir('Results_OnExcels')
end

T = 10;
h = 0.01;
xtk = 0:h:T;
n = length(xtk);
maxTimes = 5;
theta_init = 0;

alpha=0.1;
beta=0.2;
c2=1;
l=0;

%最优控制
for c1=700:50:700
    disp(c1);
[ EP,OJ,utheta,Lxt] = WAS_OptCont_Euler(T,h,maxTimes,A,theta_init,alpha,beta,c1,c2);
end
c1=700;
%静态控制
for i=1:1:10
    l=l+1;
    [ staEP,staOJ(l)] = WAS_StaCont(T,h,A,alpha,beta,c1,c2);
end
staOJ=[nan staOJ];
%%
figure;
subplot(1,2,1)
plot(xtk,utheta,'g','linewidth',2);
xlabel(c1);
ylabel('Control');

subplot(1,2,2)
stem(OJ);
hold on;
stem(staOJ);
ylabel('J(\bf{w})');

% saveas(gcf,"Fig_OnNet\sw.jpg");
%   
  
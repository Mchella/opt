clear
clc


maxTimes = 5;
theta_init = 0;

alpha=0.1;
beta=0.2;
c1=700;
c2=1;
l=1;

Ts = 1:1:10;
t=datestr(now,'yyyymmddTHHMMSS');
Igor{1,1}='alpha';
Igor{1,2}=alpha;
Igor{1,3}='c1';
Igor{1,4}=c1;
Igor{1,5}='beta';
Igor{1,6}=beta;
%最优控制
for i=1:length(Ts)
    l=l+1;
    T=Ts(i);
    disp(T);
    h = 0.01;
    xtk = 0:h:T;
    n = length(xtk);
    
    %小世界网络
    A = load('SW.txt');
    [ EP,OJ(l-1),utheta,OJS]= WAS_OptCont_Euler(T,h,maxTimes,A,theta_init,alpha,beta,c1,c2);
    Igor{l,1}=T;
    for x=1:1:maxTimes
        Igor{l,1+x}=OJS(x);
    end
    filename = ['test\test_Element_T_SW',t,'.xls'];
    xlswrite(filename,Igor);
end
%%
figure;
% subplot(1,2,1)
xlabel('T');
plot(Ts,OJ,'b','linewidth',2);
ylabel('OJ');
saveas(gcf,['Fig_OnNet\T_SW',t,'.jpg']);
%
%%
Igor = {length(Ts)+1,4};
Igor{1,1} = 'T';
Igor{1,2} = 'SW_OJ_T';
Igor{1,3}='alpha';
Igor{1,4}=alpha;
Igor{1,5}='c1';
Igor{1,6}=c1;
Igor{1,7}='beta';
Igor{1,8}=beta;
for i=2:length(Ts)+1
    Igor{i,1} = Ts(i-1);
    Igor{i,2} = OJ(i-1);
end
filename = ['Results_OnExcels\T_SW',t,'.xls'];
xlswrite(filename,Igor);

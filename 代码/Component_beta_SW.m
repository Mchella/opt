clear
clc
T = 10;
h = 0.01;
xtk = 0:h:T;
n = length(xtk);
maxTimes = 5;
theta_init = 0;

c1=1000;
c2=1;

for alpha=0.03:0.01:0.03
    
    l=1;
    betas = 0.1:0.1:1;
    t=datestr(now,'yyyymmddTHHMMSS');
    Igor{1,1}='alpha';
    Igor{1,2}=alpha;
    Igor{1,3}='c1';
    Igor{1,4}=c1;
    %最优控制
    for i=1:length(betas)
        l=l+1;
        beta = betas(i);
        disp(beta);
        %小世界网络
        A = load('SW.txt');
        [ EP,OJ(l-1),utheta,OJS]= WAS_OptCont_Euler(T,h,maxTimes,A,theta_init,alpha,beta,c1,c2);
        Igor{l,1}=beta;
        for x=1:1:maxTimes
            Igor{l,1+x}=OJS(x);
        end
        filename = ['test\test_Element_Beta_SW',t,'.xls'];
        xlswrite(filename,Igor);
    end
    %%
    figure;
    % subplot(1,2,1)
    
    plot(betas,OJ,'b','linewidth',2);
    ylabel('OJ');
    xlabel('beta');
    saveas(gcf,['Fig_OnNet\beta_SW',t,'.jpg']);
    %
    %%
    Igor = {length(betas)+1,4};
    Igor{1,1} = 'beta';
    Igor{1,2} = 'SW_OJ_beta';
    Igor{1,3}='alpha';
    Igor{1,4}=alpha;
    Igor{1,5}='c1';
    Igor{1,6}=c1;
    for i=2:length(betas)+1
        Igor{i,1} = betas(i-1);
        Igor{i,2} = OJ(i-1);
    end
    filename = ['Results_OnExcels\beta_SW',t,'.xls'];
    xlswrite(filename,Igor);
end


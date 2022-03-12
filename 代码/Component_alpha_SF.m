clear
clc
T = 10;
h = 0.01;
xtk = 0:h:T;
n = length(xtk);
maxTimes = 5;
theta_init = 0;
% 
% beta=0.1;
% c1=4000;
c1=1000;
c2=1;

for beta=0.14:0.01:0.14
    l=1;
    
    alphas = 0.1:0.1:1;
    t=datestr(now,'yyyymmddTHHMMSS');
    Igor{1,1}='beta';
    Igor{1,2}=beta;
    Igor{1,3}='c1';
    Igor{1,4}=c1;
    for i=1:length(alphas)
        l=l+1;
        alpha = alphas(i);
        disp(alpha);
        %ÎÞ±ê¶ÈÍøÂç
        load('SF20E165.mat');
        [ EP,OJ(l-1),utheta,OJS] = WAS_OptCont_Euler(T,h,maxTimes,A,theta_init,alpha,beta,c1,c2);
        Igor{l,1}=alpha;
        for x=1:1:maxTimes
            Igor{l,1+x}=OJS(x);
        end
        filename = ['test\test_Element_Alpha_SF',t,'.xls'];
        xlswrite(filename,Igor);
    end
    
    
    figure;
    % subplot(1,2,1)
    xlabel('alpha');
    plot(alphas,OJ,'b','linewidth',2);
    ylabel('OJ');
    saveas(gcf,['Fig_OnNet\alpha_SF',t,'.jpg']);
    %%
Igor = {length(alphas)+1,4};
Igor{1,1} = 'alpha';
Igor{1,2} = 'SF_OJ_alpha';
Igor{1,3}='beta';
Igor{1,4}=beta;
Igor{1,5}='c1';
Igor{1,6}=c1;
for i=2:length(alphas)+1
    Igor{i,1} = alphas(i-1);
    Igor{i,2} = OJ(i-1);
end
filename = ['Results_OnExcels\alpha_SF',t,'.xls'];
xlswrite(filename,Igor);

end


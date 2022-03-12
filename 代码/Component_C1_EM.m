clear
clc
T = 10;
h = 0.01;
xtk = 0:h:T;
n = length(xtk);
maxTimes = 5;
theta_init = 0;


beta=0.1;
c2=1;

for alpha=0.1:0.01:0.1
    for beta=0.1:0.01:0.1
        l=1;
        c1s = 1000:100:2000;
        t=datestr(now,'yyyymmddTHHMMSS');
        Igor{1,1}='alpha';
        Igor{1,2}=alpha;
        Igor{1,3}='beta';
        Igor{1,4}=beta;
        %最优控制
        for i=1:length(c1s)
            l=l+1;
            c1= c1s(i);
            disp([alpha,beta,c1]);
            %邮件网络
            [A,degree] = Get_arenas_email_Network();
            [ EP,OJ(l-1),utheta,OJS]= WAS_OptCont_Euler(T,h,maxTimes,A,theta_init,alpha,beta,c1,c2);
            Igor{l,1}=c1;
            for x=1:1:maxTimes
                Igor{l,1+x}=OJS(x);
            end
            filename = ['test\test_Element_C1_EM',t,'.xls'];
            xlswrite(filename,Igor);
        end
        %%
        figure;
        % subplot(1,2,1)
        
        plot(c1s,OJ,'b','linewidth',2);
        xlabel('c1_EM');
        ylabel('OJ');
        saveas(gcf,['Fig_OnNet\C1_EM',t,'.jpg']);
        %
        % %%
        Igor = {length(c1s)+1,4};
        Igor{1,1} = 'c1';
        Igor{1,2} = 'EM_OJ_C1';
        Igor{1,3}='alpha';
        Igor{1,4}=alpha;
        Igor{1,5}='beta';
        Igor{1,6}=beta;
        for i=2:length(c1s)+1
            Igor{i,1} = c1s(i-1);
            Igor{i,2} = OJ(i-1);
        end
        filename = ['Results_OnExcels\C1_EM',t,'.xls'];
        xlswrite(filename,Igor);
    end
end

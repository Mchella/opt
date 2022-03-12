clear
clc

T = 3;
h = 0.01;
xtk = 0:h:T;
n = length(xtk); %n=1001
maxTimes = 5;      %优化次数

x_init = 1;%控制变量：知识注入初始速率
x_Min = 0;
x_Max = 1;

Mu=0.4;
Alpha=0.2;
S=0.98*ones(n,1);
C=0.02*ones(n,1);
A=[S C];

l=0;
for lambda_init = 0.2:0.1:0.2
            w = 0.01:0.01:0.1;
            Swiki_OJ = zeros(length(w),1);
            for i=1:length(w)
                j = w(i);
                [~,~,~,AJ ] = SIPS_OptCont_Euler(T,h,maxTimes,A,x_init,x_Min,x_Max,Mu,Alpha,j);
                Swiki_OJ(i,1) = AJ;
            end
             figure('color',[1 1 1]);
            plot(T,Swiki_OJ,'r','linewidth',2);
            xlabel('\beta');
            ylabel('J^*');
            legend('wiki','location','best'); 
end
xlswrite('Results_Paras\w\w.xlsx', Swiki_OJ);
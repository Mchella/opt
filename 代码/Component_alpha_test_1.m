clear
clc
T = 10;
h = 0.01;
xtk = 0:h:T;
n = length(xtk);
maxTimes = 5;
theta_init = 0;

beta=0.1;
% c1=2500;
c1=5000;
c2=1;
l=0;

alphas = 0.1:0.1:1;
SFOJ = zeros(length(alphas),1);
SWOJ = zeros(length(alphas),1);
EMOJ = zeros(length(alphas),1);
for i=1:length(alphas)
    
    alpha = alphas(i);
    disp(i+':'+alpha);
    %无标度网络
    load('SF20E165.mat');
    [ ~,SFOJ(i),~] = WAS_OptCont_Euler(T,h,maxTimes,A,theta_init,alpha,beta,c1,c2);
    
    %小世界网络
    A = load('SW.txt');
    [ ~,SWOJ(i),~] = WAS_OptCont_Euler(T,h,maxTimes,A,theta_init,alpha,beta,c1,c2);
    
    %邮件网络
    [A,degree] = Get_arenas_email_Network();
    [ ~,EMOJ(i),~] = WAS_OptCont_Euler(T,h,maxTimes,A,theta_init,alpha,beta,c1,c2);
end


figure;
% subplot(1,2,1)
xlabel('alpha');
plot(alphas,SFOJ,'r','linewidth',2);
hold on;
plot(alphas,SWOJ,'g','linewidth',2);
hold on;
plot(alphas,EMOJ,'b','linewidth',2);
hold on;
ylabel('OJ');
legend('Scale-free','Small-world','Email','location','best');
saveas(gcf,"Fig_OnNet\alpha.jpg");
%%
Igor = {length(alphas)+1,4};
            Igor{1,1} = 'alpha';
            Igor{1,2} = 'SF_OJ_alpha';
            Igor{1,3} = 'SW_OJ_alpha';
            Igor{1,4} = 'EM_OJ_alpha';
            for i=2:length(alphas)+1
                Igor{i,1} = alphas(i-1);
                Igor{i,2} = SFOJ(i-1);
                Igor{i,3} = SWOJ(i-1);
                Igor{i,4} = EMOJ(i-1);
            end
            filename = ['Results_OnExcels\alpha.xls'];
            xlswrite(filename,Igor);

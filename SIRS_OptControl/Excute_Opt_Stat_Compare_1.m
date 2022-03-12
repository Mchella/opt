clear
clc

T = 10;
h = 0.001;
xtk = 0:h:T; n = length(xtk); %n=1001
% pn =1; % 随机策略对个数
maxTimes = 5;      %优化次数
error = 0.001;

x0 = 0;%控制变量：知识注入初始速率

xmin = 0;
xmax = 2;
omega=0.1;
mu=10;
alpha=5;
S=zeros(n,1);
C=zeros(n,1);
searchStepLength = 0.001;

converg = [];

% 记录每个节点对应的伴随变量值
lambda1=zeros(n,1);
lambda2=zeros(n,1);

x = x0*ones(n,1);    %有1个控制变量
X = xmin:searchStepLength:xmax;

% for t=1:n-1
%     dC = (func_beta1(x(t))+func_beta2(C(t)))*S(t) - alpha*C(t);
%     C(t+1)=C(t) + h*dC;
%     dS = mu - (func_beta1(x(t))+func_beta2(C(t)))*S(t) + alpha*C(t);
%     S(t+1)=S(t)+h*dS;
% end
% B = 0;
% for t=1:n-1
%     B = B + C(t)-omega*x(t);
% end
% B = B*h;
% converg = [converg B];
%最优控制
for tt=1:maxTimes
    %反向更新lamda1
    for t=n:-1:2
        diffLam1 = ( func_beta1(x(t)) + func_beta2(C(t)) )*(lambda1(t)-lambda2(t));
        lambda1(t-1) = lambda1(t)-h*diffLam1;
        diffLam2 = -1 - (alpha-func_beta2diff(C(t))*S(t))*(lambda1(t)-lambda2(t));
        lambda2(t-1) = lambda2(t)-h*diffLam2;
    end
    %更新控制变量x
    for t=1:n
        xopt = 0; optPf = -99999999999;
        for iter=1:length(X)
            pf = (lambda2(t)-lambda1(t))*S(t)*func_beta1(X(iter))-omega*X(iter);
            if pf>optPf
                optPf = pf;
                xopt = X(iter);
            end
        end
        x(t) = xopt;
    end
    % state
    for t=1:n-1
        dC = (func_beta1(x(t))+func_beta2(C(t)))*S(t) - alpha*C(t);
        C(t+1)=C(t) + h*dC;
        dS = mu - (func_beta1(x(t))+func_beta2(C(t)))*S(t) + alpha*C(t);
        S(t+1)=S(t)+h*dS;
    end
    B = 0;
    for t=1:n-1
        B = B + C(t)-omega*x(t);
    end
    B = B*h;
    converg = [converg B];
   
end

plot(x)

%静态控制
S_statics=zeros(n,1);
C_statics=zeros(n,1);
converg_statics = [];
% 下面生成初始策略
pn=100;
for i=1:pn
x_statics(:,1) = xmax*rand(n,1);
  for tt=1:maxTimes
    % state
    for t=1:n-1
        dC_statics = (func_beta1(x_statics(t))+func_beta2(C_statics(t)))*S_statics(t) - alpha*C_statics(t);
        C_statics(t+1)=C_statics(t) + h*dC_statics;
        dS_statics = mu - (func_beta1(x_statics(t))+func_beta2(C_statics(t)))*S_statics(t) + alpha*C_statics(t);
        S_statics(t+1)=S_statics(t)+h*dS_statics;
    end
    B_statics = 0;
    for t=1:n-1
        B_statics = B_statics + C_statics(t)-omega*x_statics(t);
    end
    B_statics = B_statics*h;
    %converg_statics = [converg_statics B_statics];
  end
  converg_statics=[converg_statics B_statics]
end
aj=[B converg_statics]';
xlswrite('results\example1.xlsx', x);
xlswrite('results\experiment1.xlsx', aj);
%下面求累积有效性



%最优控制
% [ ~,~,ugam,~,~,AJ ] = SIPS_OptCont_Euler(T,h,maxTimes,A,x_init,x_Min,x_Max,mu,alpha,w);
% plot(ugam)

% for p = 1:pn
%     disp(['  正在生成第',num2str(p),'个随机策略对...']);
%     % 下面生成初始策略
%     x_statics(:,1) = x_Max*rand(n,1);
%     %静态控制
%     [ ugam1,AL1,AC1,AJ1 ] = SIPS_StaCont(T,h,A,x_statics,mu,alpha,w);
%     A_AJ1(:,p)=AJ1;
% end

% cp=[AJ A_AJ1];
% xlswrite('Results_OnWiki\Opt_stat_compare\Excels\x.xlsx', ugam);
%xlswrite('Results_OnWiki\Opt_stat_compare\Excels\compare.xlsx', cp);



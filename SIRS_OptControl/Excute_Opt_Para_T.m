clear
clc
BB=[];
hjm=2:2:40;
for i=1:length(hjm)
    h = 0.0001;
% pn =1; % 随机策略对个数
maxTimes = 5;      %优化次数
error = 0.001;
x0 = 0;%控制变量：知识注入初始速率
xmin = 0;
xmax = 2;
omega=0.1;
mu=10;
alpha=5;
searchStepLength = 0.001;
converg = [];
X = xmin:searchStepLength:xmax;
%最优控制

    TT= hjm(i);
    xtk = 0:h:TT; n = length(xtk); %n=1001
    S=zeros(n,1);
    C=zeros(n,1);
    % 记录每个节点对应的伴随变量值
    lambda1=zeros(n,1);
    lambda2=zeros(n,1);
    x = x0*ones(n,1);    %有1个控制变量
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
    
end
BB(i)=B;
end
xlswrite('results\experiment4.xlsx', BB');

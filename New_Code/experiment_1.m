function eReturn=experiment_1()
clc
clear
T=5;%时间
S_0=[100,100,0];

theta_MIN=0;
theta_MAX=5;
mu=1;
delta=0.01;
alpha=0.05;
%beta_2;
gamma=0.05;
omega_1=3;
omega_2=5;
maxTime=40;
h=0.01;
stepNum=0:h:T;
n=length(stepNum);
%x0=theta_MAX;

searchStepLength=0.01;
lambda_1=zeros(n,1);
lambda_2=zeros(n,1);
lambda_3=zeros(n,1);
I=S_0(1)*ones(n,1);
U=S_0(2)*ones(n,1);
W=S_0(3)*ones(n,1);
%sum=zeros(n,1);


theta = theta_MAX*ones(n,1);    %有1个控制变量
Theta = theta_MIN:searchStepLength:theta_MAX;

%无控制的计算
for t=1:n-1
    
    dI=mu-alpha*I(t)*W(t)-delta*I(t);
    I(t+1)=I(t)+h*dI;
    dU=alpha*I(t)*W(t)-beta_1(theta(t))*U(t)+gamma*W(t)-delta*U(t);
    U(t+1)=U(t)+h*dU;
    dW=beta_1(theta(t))*U(t)-gamma*W(t)-delta*W(t);
    W(t+1)=W(t)+h*dW;
        if(I(t+1)<0)
            I(t+1)=0;
        end
         if(U(t+1)<0)
            U(t+1)=0;
         end
        if(W(t+1)<0)
            W(t+1)=0;
        end
end

J=0;%目标函数
for t=1:n
    dJ=omega_1*I(t)+omega_2*(U(t)+W(t));
    J=J+h*dJ;
end

eReturn=[J];

%最优控制

for tt=1:maxTime
    disp(['最优控制第',num2str(tt),'次迭代']);
    %反向计算lambda
    for t=n:-1:2
        dLambda_1=(alpha*W(t)+delta)*lambda_1(t)-alpha*W(t)*lambda_2(t)-omega_1;
        lambda_1(t-1)=lambda_1(t)-h*dLambda_1;
        dLambda_2=(beta_1(theta(t))+delta)*lambda_2(t)-beta_1(theta(t))*lambda_3(t)-omega_2;
        lambda_2(t-1)=lambda_2(t)-h*dLambda_2;
        dLambda_3=alpha*I(t)*lambda_1(t)-(alpha*I(t)+gamma)*lambda_2(t)+(gamma+delta)*lambda_3(t)+theta(t)-omega_2;
        lambda_3(t-1)=lambda_3(t)-h*dLambda_3;   
    end
        
    %更新控制变量
    %thetaOpt=0;
    %beta_1=[];
    for t=1:n
        optPf = -99999999999;
        for iter=1:length(Theta)
            
            pf = (lambda_3(t)-lambda_2(t))*U(t)*beta_1(Theta(iter))-W(t)*Theta(iter);
            if pf>=optPf
                optPf = pf;
                theta(t) = Theta(iter);
            end
        end
        %theta(t) = thetaOpt; 
    end
    %theta_1=[theta_1 theta];
    
    %更新状态
    for t=1:n-1
    
        dI=mu-alpha*I(t)*W(t)-delta*I(t);
        I(t+1)=I(t)+h*dI;
        dU=alpha*I(t)*W(t)-beta_1(theta(t))*U(t)+gamma*W(t)-delta*U(t);
        U(t+1)=U(t)+h*dU;
        dW=beta_1(theta(t))*U(t)-gamma*W(t)-delta*W(t);
        W(t+1)=W(t)+h*dW;
        if(I(t+1)<0)
            I(t+1)=0;
        end
         if(U(t+1)<0)
            U(t+1)=0;
         end
         if(W(t+1)<0)
            W(t+1)=0;
        end
        
    
    end
    J=0;%目标函数
    for t=1:n
        dJ=omega_1*I(t)+omega_2*(U(t)+W(t));
        J=J+h*dJ;
    end
    eReturn=[eReturn J];
    
      
end

eRandReturn=[];
I_statics=S_0(1)*ones(n,1);
U_statics=S_0(2)*ones(n,1);
W_statics=S_0(3)*ones(n,1);

for nn=1:100
    
theta_statics=zeros(n,1);
theta_statics(:,1)=theta_MAX*rand(n,1);
for t=1:n
    
        dI_statics=mu-alpha*I_statics(t)*W_statics(t)-delta*I_statics(t);
        I_statics(t+1)=I_statics(t)+h*dI_statics;
        dU_statics=alpha*I_statics(t)*W_statics(t)-beta_1(theta_statics(t))*U_statics(t)+gamma*W_statics(t)-delta*U_statics(t);
        U_statics(t+1)=U_statics(t)+h*dU_statics;
        dW_statics=beta_1(theta_statics(t))*U_statics(t)-gamma*W_statics(t)-delta*W_statics(t);
        W_statics(t+1)=W_statics(t)+h*dW_statics;
        if(I_statics(t+1)<0)
            I_statics(t+1)=0;
        end
         if(U_statics(t+1)<0)
            U_statics(t+1)=0;
         end
         if(W_statics(t+1)<0)
            W_statics(t+1)=0;
         end
end
    J_R=0;%目标函数
    for t=1:n
        dJ_R=omega_1*I_statics(t)+omega_2*(U_statics(t)+W_statics(t));
        J_R=J_R+h*dJ_R;
    end
    eRandReturn=[eRandReturn J_R];
end
sand=4000*ones(1,5);
    Result=[eReturn(maxTime)  eRandReturn];

% 
% 
% subplot(2,3,1)
%     plot(theta)
%     subplot(2,3,2)
%     plot(I)
%     subplot(2,3,3)
%     plot(U)
%     subplot(2,3,4)
%     plot(W)
%     subplot(2,3,5)
%     plot(eReturn)
%     subplot(2,3,6)
%     plot(Result)
%    xlswrite('results\experiment1.xlsx', Result);


minIn=min(Result);
maxIn=max(eRandReturn);





plot(10,Result(1),'r.','MarkerSize',20);
for i=1:100
    hold on
    plot(i+20,eRandReturn(i),'b.','MarkerSize',20)
end


for i=1:100
    %hold on
    line([i+20,i+20],[minIn-100,eRandReturn(i)])
end

line([10,10],[minIn-100,Result(1)],'Color','red')

xlabel('$\theta$','interpreter','latex');
ylabel('$J(\theta)$','interpreter','latex');

set(gca, 'xTick', [0:10:120]);

set(gca,'XTickLabel',{'','\theta_1^{*}','\theta_1','','','','…','','','','','','\theta_{100}'})
set(gca,'color','w')    
%set(gca, 'FontSize', 18);

line([0,120],[minIn,minIn],'Color','black','LineStyle','--');

line([0,120],[maxIn,maxIn],'Color','black','LineStyle','--');
axis([0 120,minIn-100,max(Result)+100]);

% xlim([0, 125])
% ylim([minIn-100, max(Result)+100])
%title('(b)')
set(gca,'FontSize',19);
     %set(gca,'FontSize',50);
    print fig\experiment_1.eps -depsc2 -r600

end



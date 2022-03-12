function eReturn=ex_3()
clc
clear
 
T=5;%时间
S_0=[100,100,0];

theta_MIN=0;
theta_MAX=5;
mu=1;
delta=0.1;
alpha=0.01;
%beta;
gamma=0.05;
omega_1=20;
omega_2=50;
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
% for t=1:n-1
%     
%     dI=mu-alpha*I(t)*W(t)-delta*I(t);
%     I(t+1)=I(t)+h*dI;
%     dU=alpha*I(t)*W(t)-beta(theta(t))*U(t)+gamma*W(t)-delta*U(t);
%     U(t+1)=U(t)+h*dU;
%     dW=beta(theta(t))*U(t)-gamma*W(t)-delta*W(t);
%     W(t+1)=W(t)+h*dW;
%         if(I(t+1)<0)
%             I(t+1)=0;
%         end
%          if(U(t+1)<0)
%             U(t+1)=0;
%          end
%         if(W(t+1)<0)
%             W(t+1)=0;
%         end
% end
% 
% J=0;%目标函数
% for t=1:n
%     dJ=omega_1*I(t)+omega_2*(U(t)+W(t));
%     J=J+h*dJ;
% end
% 
eReturn=[];

%最优控制

PF=[];
for tt=1:maxTime
    disp(['最优控制第',num2str(tt),'次迭代']);
    %反向计算lambda
    for t=n:-1:2
        dLambda_1=(alpha*W(t)+delta)*lambda_1(t)-alpha*W(t)*lambda_2(t)-omega_1;
        lambda_1(t-1)=lambda_1(t)-h*dLambda_1;
        dLambda_2=(beta_3(theta(t))+delta)*lambda_2(t)-beta_3(theta(t))*lambda_3(t)-omega_2;
        lambda_2(t-1)=lambda_2(t)-h*dLambda_2;
        dLambda_3=alpha*I(t)*lambda_1(t)-(alpha*I(t)+gamma)*lambda_2(t)+(gamma+delta)*lambda_3(t)+theta(t)-omega_2;
        lambda_3(t-1)=lambda_3(t)-h*dLambda_3;   
    end
        
    %更新控制变量
    %thetaOpt=0;
    
    for t=1:n
        optPf = -inf;   
        for iter=1:length(Theta)
            pf = (lambda_3(t)-lambda_2(t))*U(t)*beta_3(Theta(iter))-W(t)*Theta(iter);
            PF(tt,t)=pf;
            if pf>optPf
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
        
        dU=alpha*I(t)*W(t)-beta_3(theta(t))*U(t)+gamma*W(t)-delta*U(t);
        U(t+1)=U(t)+h*dU;
        dW=beta_3(theta(t))*U(t)-gamma*W(t)-delta*W(t);
        W(t+1)=W(t)+h*dW;
%         if(I(t+1)<0)
%             I(t+1)=0;
%         end
%          if(U(t+1)<0)
%             U(t+1)=0;
%          end
%          if(W(t+1)<0)
%             W(t+1)=0;
%         end
        
    
    end
    J=0;%目标函数
    for t=1:n
        dJ=omega_1*I(t)+omega_2*(U(t)+W(t));
        J=J+h*dJ;
    end
    eReturn=[eReturn J];
    
      
end
%xlswrite('results\example3.xlsx', theta);
% subplot(2,3,1)
    plot(theta)
%     subplot(2,3,2)
%     plot(I)
%     subplot(2,3,3)
%     plot(U)
%     subplot(2,3,4)
%     plot(W)
%     subplot(2,3,5)
%     plot(eReturn)
    
    
    xlabel('$t$','interpreter','latex');
    axis([0 500,-1,6]);
    ylabel('$\theta_3^*$','interpreter','latex');
    %set(gca, 'xTick', [0:100:500]);
    %set(gca, 'yTick', [0:1:7]);
   % set(gca, 'xTick', [0:10:120]);
    %set(gca,'FontSize',30);

    set(gca,'XTickLabel',{'0','1','2','3','4','5'})

     set(gca,'FontSize',25);
     %set(gca,'FontSize',55);
    print fig\ex_3.eps -depsc2 -r600




end



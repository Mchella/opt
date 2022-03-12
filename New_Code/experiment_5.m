function eNumReturn=experiment_5()
clc
clear

eNumReturn=[];
for num=1:10   
 T=5;%ʱ�� 
   
    disp(['theta=',num2str(num),'�ε���']);
S_0=[100,100,0];

theta_MIN=0;
theta_MAX=num;
mu=2;
delta=0.01;
alpha=0.05;
%beta_e_4;
gamma=0.05;
omega_1=3;
omega_2=5;
maxTime=15;
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


theta = theta_MAX*ones(n,1);    %��1�����Ʊ���
Theta = theta_MIN:searchStepLength:theta_MAX;

%�޿��Ƶļ���
for t=1:n-1
    
    dI=mu-alpha*I(t)*W(t)-delta*I(t);
    I(t+1)=I(t)+h*dI;
    dU=alpha*I(t)*W(t)-beta_e_4(theta(t))*U(t)+gamma*W(t)-delta*U(t);
    U(t+1)=U(t)+h*dU;
    dW=beta_e_4(theta(t))*U(t)-gamma*W(t)-delta*W(t);
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

J=0;%Ŀ�꺯��
for t=1:n
    dJ=omega_1*I(t)+omega_2*(U(t)+W(t));
    J=J+h*dJ;
end

eReturn=[J];

%���ſ���

for tt=1:maxTime
    disp(['���ſ��Ƶ�',num2str(tt),'�ε���']);
    %�������lambda
    for t=n:-1:2
        dLambda_1=(alpha*W(t)+delta)*lambda_1(t)-alpha*W(t)*lambda_2(t)-omega_1;
        lambda_1(t-1)=lambda_1(t)-h*dLambda_1;
        dLambda_2=(beta_e_4(theta(t))+delta)*lambda_2(t)-beta_e_4(theta(t))*lambda_3(t)-omega_2;
        lambda_2(t-1)=lambda_2(t)-h*dLambda_2;
        dLambda_3=alpha*I(t)*lambda_1(t)-(alpha*I(t)+gamma)*lambda_2(t)+(gamma+delta)*lambda_3(t)+theta(t)-omega_2;
        lambda_3(t-1)=lambda_3(t)-h*dLambda_3;   
    end
        
    %���¿��Ʊ���
    %thetaOpt=0;
    %beta_e_4=[];
    for t=1:n
        optPf = -99999999999;
        for iter=1:length(Theta)
            
            pf = (lambda_3(t)-lambda_2(t))*U(t)*beta_e_4(Theta(iter))-W(t)*Theta(iter);
            if pf>=optPf
                optPf = pf;
                theta(t) = Theta(iter);
            end
        end
        %theta(t) = thetaOpt; 
    end
    %theta_1=[theta_1 theta];
    
    %����״̬
    for t=1:n-1
    
        dI=mu-alpha*I(t)*W(t)-delta*I(t);
        I(t+1)=I(t)+h*dI;
        dU=alpha*I(t)*W(t)-beta_e_4(theta(t))*U(t)+gamma*W(t)-delta*U(t);
        U(t+1)=U(t)+h*dU;
        dW=beta_e_4(theta(t))*U(t)-gamma*W(t)-delta*W(t);
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
    J=0;%Ŀ�꺯��
    for t=1:n
        dJ=omega_1*I(t)+omega_2*(U(t)+W(t));
        J=J+h*dJ;
    end
    eReturn=[eReturn J];
    
      
end
   eNumReturn=[eNumReturn J];
end
xlswrite('results\experiment5.xlsx', eNumReturn);
% subplot(2,3,1)
    plot([1:1:10],eNumReturn,'ro-')
    set(gcf,'color','w')  
    xlabel('','Interpreter','latex','String','$\bar \theta$')
    %xlabel('$\overline \theta$');
    ylabel('','Interpreter','latex','String','$J(\theta^{*}_{\bar \theta})$')
    set(gca,'FontSize',19);
    %set(gca,'FontSize',20);
    print fig\experiment_5.eps -depsc2 -r600
%    subplot(1,3,2)
%    plot(I)
%    subplot(2,3,3)
%    plot(U)
%     subplot(2,3,4)
%    plot(W)
%    subplot(2,3,5)
%    plot(eReturn)
%     





end



function [ EP,OJ,utheta,OJS] = WAS_OptCont_Euler(T,h,maxTimes,A,theta_init,alpha,beta,c1,c2)
global L M  Active Silent theta

[N,~] = size(A);

% h:����,dt
% n:ʱ�䳤��
% N:�ڵ����

t0 = 0:h:T;
n = length(t0);

% %��ʼ��E0
Active =zeros(n,N);
% Silent= zeros(n,N);
% Active =0.05*ones(n,N);
Silent= 0.1*ones(n,N);

L = zeros(n,N);         %lamuda
M = zeros(n,N);         %miu
L(n,:) = 0;
M(n,:) = 0;
theta =theta_init*ones(n,1);    %��ʼ��theta��һ��Ϊ1

OJS=zeros(maxTimes,1);

flag=true;
for tt=1:maxTimes
    if(flag==false)
        break;
    end
    disp(['���ſ��Ƶ�',num2str(tt),'�ε���']);
    for t=1:n-1
        %����Active
        for i=1:N
            y = Active (t,i);
            result = y + h*fA(t,i,alpha);
            if result>1
                disp('A>1');
                Active (t+1,i)=1;
            elseif result<0
                disp('A<0');
                Active (t+1,i)=0;
            else
                Active (t+1,i)=result;
            end
        end
        %����Silent
        for i=1:N
            y = Silent(t,i);
            result = y + h*fS(t,i,A,alpha,beta);
            if result>1
                disp('S>1');
                Silent(t+1,i)=1;
            elseif result<0
                disp('S<0');
                Silent(t+1,i)=0;
            else
                Silent(t+1,i)=result;
            end
        end
    end
    %����lambda��mu
    for t=n:-1:2
        % lambda
        for i=1:N
            y = L(t,i);
            L(t-1,i) = y-h*fL(t,i,A,alpha,beta,c1,c2);
        end
        % miu
        for i=1:N
            y = M(t,i);
            M(t-1,i) = y-h*fM(t,i,A,beta,c1,c2);
        end
    end
    
    %����theta
    for t=2:n
        % ������һ�����⣬sum1�����Ƿ���ѭ���⻹��ѭ����
        sum1=0;
        sum2=0;
        sum3=0;
        sum4=0;
        for i=1:N
            for j=1:N
                sum1=sum1+(1-Silent(t,i)-Active(t,i))*A(i,j)*Active(t,j);
                sum2=sum2+Silent(t,i)+Active(t,i);
                sum3=sum3+A(i,j)*Active(t,j)*L(t,i)*(1-Silent(t,i)-Active(t,i));
            end
        end
        sum4=sum4+c1*beta*sum1-c2*sum2+beta*sum3;
        if sum4<0
            theta(t) = 0;
        else
            theta(t) = 1;
        end
    end
    
    EP = zeros(n,1);
    Lxt = zeros(n,1);
    
    for t=1:n
        sum1=0;
        for i=1:N
            for j=1:N
                sum1=sum1+(1-Active(t,i)-Silent(t,i))*A(i,j)*Active(t,j);
            end
        end
        sum2=0;
        for i=1:N
            sum2=sum2+Active(t,i)+Silent(t,i);
        end
        Lxt(t) =c1*beta*theta(t)*sum1-c2*sum2*theta(t);
    end
    EP(1) = Lxt(1);
    for j=2:n
        EP(j) = EP(j-1)+Lxt(j);
    end
    OJS(tt,1) = EP(n);
    if(EP(n)==0)
        flag=false;
    end    
end
disp(OJS);

EP = zeros(n,1);
Lxt = zeros(n,1);

% if(flag==false)
%     OJ=EP(n);
%     utheta=theta;
%     return;
% end

for t=1:n
    sum1=0;
    for i=1:N
        for j=1:N
            sum1=sum1+(1-Active(t,i)-Silent(t,i))*A(i,j)*Active(t,j);
        end
    end
    sum2=0;
    for i=1:N
        sum2=sum2+Active(t,i)+Silent(t,i);
    end
    Lxt(t) =c1*beta*theta(t)*sum1-c2*sum2*theta(t);
end
EP(1) = Lxt(1);
for j=2:n
    EP(j) = EP(j-1)+Lxt(j);
end
OJ=EP(n);
utheta=theta;
end





function [ EP,OJ] = WAS_StaCont(T,h,A,alpha,beta,c1,c2)

global  Active_static Silent_static theta_static;
[N,~] = size(A);
% h:����,dt
% n:ʱ�䳤��
% N:�ڵ����
 
t0 = 0:h:T;    
n = length(t0); 

% ��ʼ��E0
Active_static = zeros(n,N);
% Silent_static = zeros(n,N);
% Active_static =0.1*ones(n,N);
Silent_static = 0.1*ones(n,N);

% theta_static = static_theta*ones(n,1);

% ��ʼ��theta_static,�������
theta_static =rand(n,1);


disp('���ڴ���̬����');
for t=1:n-1
    %����Active
    for i=1:N
        y = Active_static(t,i);
        result = y + h*fA_static(t,i,alpha);
        if result>1
            disp('A>1');
            Active_static(t+1,i)=1;
        elseif result<0
            disp('A<0');
            Active_static(t+1,i)=0;
        else
            Active_static(t+1,i)=result;
        end
    end
    %����Slient
    for i=1:N
        y = Silent_static(t,i);
        result = y + h*fS_static(t,i,A,alpha,beta);
        if result>1
            disp('S>1');
            Silent_static(t+1,i)=1;
        elseif result<0
            disp('S<0');
            Silent_static(t+1,i)=0;
        else
            Silent_static(t+1,i)=result;
        end
    end
end


%����������
EP = zeros(n,1);
Lxt = zeros(n,1);
for t=1:n
    sum1=0;
    for i=1:N
        for j=1:N
            sum1=sum1+(1-Active_static(t,i)-Silent_static(t,i))*A(i,j)*Active_static(t,j);
        end
    end
    sum2=0;
    for i=1:N
        sum2=sum2+Active_static(t,i)+Silent_static(t,i);
    end
    Lxt(t) =c1*beta*theta_static(t)*sum1-c2*sum2*theta_static(t);
end
    EP(1) = Lxt(1);
    for j=2:n
        EP(j) = EP(j-1)+Lxt(j);
    end
    OJ=EP(n);
end
 
 
 
 
 
 



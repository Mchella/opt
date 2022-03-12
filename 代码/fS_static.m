function [ dS ] = fS_static(t,i,A,alpha,beta)
global Silent_static Active_static theta_static
[N,~]=size(A);
sum1=0;
for j=1:N
    sum1=sum1+A(i,j)*Active_static(t,j);
end

dS=beta*theta_static(t)*(1-Silent_static(t,i)-Active_static(t,i))*sum1-alpha*Silent_static(t,i);
end


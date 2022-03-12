function [ dS ] = fS(t,i,A,alpha,beta)

global Silent Active theta
[N,~]=size(A);
sum1=0;
for j=1:N
    sum1=sum1+A(i,j)*Active(t,j);
end

dS=beta*theta(t)*(1-Silent(t,i)-Active(t,i))*sum1-alpha*Silent(t,i);
    
end


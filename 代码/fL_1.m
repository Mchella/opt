function [ dL ] = fL(t,i,A,alpha,beta,c1,c2)

global Active theta L M
 [N,~]=size(A);
sum1=0;
for j=1:N
    sum1=sum1+A(i,j)*Active(t,j);
end

dL=c1*beta*theta(t)*sum1+c2*theta(t)+(beta*theta(t)*sum1+alpha)*L(t,i)-alpha*M(t,i);

end
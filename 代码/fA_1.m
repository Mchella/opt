function [ dA ] = fA(t,i,alpha)
% ����
% t -- ����ʱ��
% i -- �ڵ���
% A -- �ڽӾ���

global   Silent;
% [N,~]=size(A);
% sum1=0;
% for j = 1:N
%     sum1 = sum1+A(i,j)*Positive(t,j);
% end
% dA=alpha*sum1*(1-Adopting(t,i)-Positive(t,i))-gamma1*Adopting(t,i)-(beta1+beta2*udelta(t))*Adopting(t,i);
dA=alpha*Silent(t,i);
end
function [ dA ] = fA(t,i,alpha)
% 参数
% t -- 迭代时刻
% i -- 节点编号
% A -- 邻接矩阵

global   Silent;
% [N,~]=size(A);
% sum1=0;
% for j = 1:N
%     sum1 = sum1+A(i,j)*Positive(t,j);
% end
% dA=alpha*sum1*(1-Adopting(t,i)-Positive(t,i))-gamma1*Adopting(t,i)-(beta1+beta2*udelta(t))*Adopting(t,i);
dA=alpha*Silent(t,i);
end
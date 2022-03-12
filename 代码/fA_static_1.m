function [ dA ] = fA_static(t,i,alpha)
% 参数
% t -- 迭代时刻
% i -- 节点编号
% di -- 归一化的节点度


global   Silent_static;
% [N,~]=size(A);
% sum1=0;
% for j = 1:N
%    sum1 = sum1+A(i,j)*Positive_static(t,j);
% end
dA=alpha*Silent_static(t,i);
end

function [ dA ] = fA_static(t,i,alpha)
% ����
% t -- ����ʱ��
% i -- �ڵ���
% di -- ��һ���Ľڵ��


global   Silent_static;
% [N,~]=size(A);
% sum1=0;
% for j = 1:N
%    sum1 = sum1+A(i,j)*Positive_static(t,j);
% end
dA=alpha*Silent_static(t,i);
end

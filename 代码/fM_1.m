function [ dM ] = fM(t,i,A,beta,c1,c2)
% 参数
% t -- 迭代时刻
% i -- 节点编号
% A -- 邻接矩阵

global   Active Silent L theta;
[N,~]=size(A);
    sum1=0;
    sum2=0;
    sum3=0;
    for j=1:N
        sum1=sum1+A(i,j)*(2*Active(t,j)+Silent(t,j)-1);
        % 注意，这里到底是 A(j,i)还是 A(i,j)
        sum2=sum2+A(i,j)*Active(t,j);
        sum3=sum3+L(t,j)*A(i,j)*(1-Active(t,j)-Silent(t,j));
    end
    dM=c1*beta*theta(t)*sum1+c2*theta(t)+beta*theta(t)*(L(t,i)*sum2-sum3);
end
function [A] = smallworld(N,m,p,M)
%小世界网络生成算法  20.10.26
% N个节点，换边概率为p的小世界网络  
%   默认m=2 即在Lattice基础上随机换边  by timothy 20170322

% 生成一个规则网络 m*2正则图
edges = 0;
% A=sparse(N,N);
A = M;
if mod(m,2)~=0
    disp('m must be even');
    return;
end
for i=1:N
    for j=1:m
        source = i;
        if i+j>N
            target = mod(i+j,N);
        else
            target = i+j;
        end
        A(source,target) = 1;
        A(target,source) = 1;
        % 加入边集
        edges = edges + 1;
        edgeSet(edges).s = source; 
        edgeSet(edges).t = target;
    end
end

%开始随机换边
for i=1:edges
    tempR = rand(1);
    if tempR <= p
        while 1 % 删边
            deleteNO = randi(edges,1);
            source = edgeSet(deleteNO).s;
            target = edgeSet(deleteNO).t;
            if sum(A(source,:))>1 && sum(A(target,:))>1 % 不能删度为1的节点的边
                A(source,target) = 0;
                A(target,source) = 0;
                break;
            end     
        end
        while 1 % 加边
             source = randi(N,1);
             target = randi(N,1);
             if source~=target && A(source,target)==0 % 不能加环和平行边
                 A(source,target) = 1;
                 A(target,source) = 1;
                 edgeSet(deleteNO).s = source;
                 edgeSet(deleteNO).t = target;
                 break;
             end
        end        
    end
end

end


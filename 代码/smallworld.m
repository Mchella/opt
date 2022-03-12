function [A] = smallworld(N,m,p,M)
%С�������������㷨  20.10.26
% N���ڵ㣬���߸���Ϊp��С��������  
%   Ĭ��m=2 ����Lattice�������������  by timothy 20170322

% ����һ���������� m*2����ͼ
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
        % ����߼�
        edges = edges + 1;
        edgeSet(edges).s = source; 
        edgeSet(edges).t = target;
    end
end

%��ʼ�������
for i=1:edges
    tempR = rand(1);
    if tempR <= p
        while 1 % ɾ��
            deleteNO = randi(edges,1);
            source = edgeSet(deleteNO).s;
            target = edgeSet(deleteNO).t;
            if sum(A(source,:))>1 && sum(A(target,:))>1 % ����ɾ��Ϊ1�Ľڵ�ı�
                A(source,target) = 0;
                A(target,source) = 0;
                break;
            end     
        end
        while 1 % �ӱ�
             source = randi(N,1);
             target = randi(N,1);
             if source~=target && A(source,target)==0 % ���ܼӻ���ƽ�б�
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


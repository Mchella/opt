function  [matrix,degree]  = Get_arenas_email_Network( )
%100个节点的网络导入 20.10.26 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% load arenas-email.txt
load short_social.txt
% Edges = arenas_email;
Edges = short_social;
N = max(max(Edges));
matrix = zeros(N,N);
l = length(Edges);
for i=1:l
    edge = Edges(i,:);
    matrix(edge(1),edge(2)) = 1;
    matrix(edge(2),edge(1)) = 1;
end
degree = zeros(N,1);
for i=1:N
    degree(i) = sum(matrix(i,:));
end

end


%检验不同网络中的合适参数得到的最优控制收敛
clear
clc
A = load('SW.txt'); % 小世界网络邻接矩阵导入A
if ~exist('test')
    mkdir('test')
end

T = 10;
h = 0.01;
xtk = 0:h:T;
n = length(xtk);
maxTimes = 6;
theta_init = 0;

% alpha=0.1;
% beta=0.1;
c2=1;
x=1;
t=datestr(now,'yyyymmddTHHMMSS');
Igor{x,1}='alpha';
Igor{x,2}='beta';
Igor{x,3}='c1';
Igor{x,4}='OJS';
x=x+1;
%最优控制
for c1=500:100:700
    for alpha=0.01:0.01:0.2
        for beta=0.05:0.01:0.2
            disp([c1,alpha,beta]);
            [ EP,OJ,utheta,OJS] = WAS_OptCont_Euler(T,h,maxTimes,A,theta_init,alpha,beta,c1,c2);
            Igor{x,1}=alpha;
            Igor{x,2}=beta;
            Igor{x,3}=c1;
            for i=1:1:maxTimes
                Igor{x,3+i}=OJS(i);
            end
            if OJS(maxTimes)==0
                continue;
            end
        filename = ['test\test_final_SW',t,'.xls'];
        xlswrite(filename,Igor);
        x=x+1;
        end  
    end

end

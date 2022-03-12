function [l2,l1,ux,AL,AC,AJ ] = SIPS_OptCont_Euler(T,h,maxTimes,A,x_init,x_Min,x_Max,mu,alpha,w)

%本程序使用一阶欧拉法求解微分方程组

%参数
% T --- 迭代时间
% h --- Euler法的步长
% maxTimes --- 循环求解的次数

%返回值
% CP --- 累积有效性
% ugam -- 最优控制
% AL -- 平均损失
% AC -- 平均能源消耗
% AJ -- AL + AC

global S C x lamda1 lamda2;




function [l2,l1,ux,AL,AC,AJ ] = SIPS_OptCont_Euler(T,h,maxTimes,A,x_init,x_Min,x_Max,mu,alpha,w)

%������ʹ��һ��ŷ�������΢�ַ�����

%����
% T --- ����ʱ��
% h --- Euler���Ĳ���
% maxTimes --- ѭ�����Ĵ���

%����ֵ
% CP --- �ۻ���Ч��
% ugam -- ���ſ���
% AL -- ƽ����ʧ
% AC -- ƽ����Դ����
% AJ -- AL + AC

global S C x lamda1 lamda2;




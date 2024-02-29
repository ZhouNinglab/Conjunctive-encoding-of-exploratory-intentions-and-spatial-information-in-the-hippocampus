function [lapS] = rand_label(explap)
%RAND_LABEL �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

nlap1 = numel(explap.explap1); 
nlap2 = numel(explap.explap2);
lapnum = nlap1 + nlap2;

explap1s = randperm(lapnum,nlap1);
explap2s = 1:lapnum;
explap2s(explap1s)=[];
lapS.explap1 = explap1s;
lapS.explap2 = explap2s;


end


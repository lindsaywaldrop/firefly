function [c]=getHairSniff(n,pathbase,plotn)
% getHairSniff.m
%
% Script for analyzing hair data from Shilpa's sniffing code. 
%

nn=sprintf('%0.4d',n);

load([pathbase,'initdata_',num2str(nn),'.mat']);
load([pathbase,'hairs_c_',num2str(nn),'.mat']);

s=floor(t_final_flick/dt_flick/40)-3;

sumHairs=zeros(s,3);

for i = 1:s
 
    val=eval(['hairs_c_',num2str(i)]);
    
    eval(['sumHairs(i,1)=sum(hairs_c_',num2str(i),'{1,1});'])
    eval(['sumHairs(i,2)=sum(hairs_c_',num2str(i),'{1,2});'])
    eval(['sumHairs(i,3)=sum(hairs_c_',num2str(i),'{1,3});'])
    
end

if plotn==1
%    plot(list_print_times,sumHairs)
else
    
end

c = sum(sumHairs(s,:));

% calcHairs.m
%
% Script for extracting and summing all concentrations from Shilpa's sniffing
% code. Returns a single value for each simulation. 
%

clear all
close all

pathbase=('/Users/Bosque/Dropbox/EntIBAMR/results/completed_water/');

%list=[1,8,23,83,146,156,158,175,199,244,247,297,299,303,316,356,431,443,...
%517,533,552,557,562,601,602,608,618,625,633,640,666,685,690,695,703,717,800,...
%861,931,946,998,1050,1116];
%conc = zeros(length(list),1)
y = 1233; % number of simulations to analyze
conc = zeros(y,1);

missing=0;

for u=1:y
%for k=1:length(list)
    %u=list(k); 
    disp(['Simulation: ',num2str(u)])
    try
    conc(u)=getHairSniff(u,pathbase,0);
    disp(num2str(conc(u)))
    catch
    missing(u)=u;
    end
end


%plot(conc)

csvwrite('ConcWater2018.csv',conc)

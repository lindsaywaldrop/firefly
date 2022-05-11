% convert_entcode.m
%	This script imports SAMRAI data into Matlab, interpolates to a
%	uniform grid. Only the final time step (set by user) is imported.
%
%   Note: This script relies on Eric Tytell's Matlab code for importing samrai data! 
%	This package needs to be installed for this to work properly!
% 
% 	GridSize = IBAMR fluid grid size
%	final_time = final time step of simulation, or time step of interest.
% 	n = total number of simulations
%	pathbase = where runs are located

close all
clear all
clc

GridSize = 4096;        %Size of the finest grid
final_time=25000;      %Final time step to be included in data analysis
n=1233; %Number of simulations

pathbase='/Volumes/HelmsDeep/IBAMR/entcode/usedruns/';

%for i=1:n
    i=1233
    clear V Vinterp x y
    file=['viz_IB2d',num2str(i)];
    path1=[pathbase,file];
    
    disp(['Simulation number: ',num2str(i)])
    
    %Make this the pathname for the visit data
	path_name1 = [path1,'/visit_dump.',num2str(final_time)];
    
    %This reads the samrai data and interpolates it.
    [x,y,Vinterp,V] = importsamrai(path_name1,'interpolaten',[GridSize GridSize]);
    
    %Shilpa's code goes here. 
    
    
    % Saves data relevant to next step in workflow.
    % Note: to tighten workflow, we won't be saving this step. 
    save([file,'.mat'],'V','Vinterp','x','y','GridSize','final_time')

%end

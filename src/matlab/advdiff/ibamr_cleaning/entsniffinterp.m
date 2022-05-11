function entsniffinterp(i,files,pathbase1,GridSize,final_time)
%
%
%

%This reads the samrai data and interpolates it.
[x,y,Vinterp,V] = importsamrai([pathbase1,'viz_IB2d',files{i},'/visit_dump.', final_time],'interpolaten',[GridSize GridSize]);

% Saves data relevant to next step in workflow.                              
save([pathbase1,'viz_IB2d',files{i},'.mat'],'V','Vinterp','x','y','GridSize','final_time')

clear x y Vinterp V

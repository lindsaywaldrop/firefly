function [hairNum] = convert_hairdata(pathbase_data, Species, run_id)
%
%

run_number = run_id

hairdata = csvread(strcat(pathbase_data, 'csv-files/',Species,...
    '/',Species,'_',num2str(run_number),'.csv'),1);
numPoints = hairdata(end,1) + hairdata(end,2);
hairNum = length(hairdata);

hair_vertices = dlmread(strcat(pathbase_data,'vertex-files/',Species,...
        '/',Species,'_', num2str(run_id),'.vertex'),' ');
hair_vertices = hair_vertices(2:numPoints+1, :);
    
for h = 1:hairNum
    disp(num2str(h))
    %disp('centers')
    %eval(['p.hair',num2str(h),'Centerx = hairdata(2,h+1);']);
    %eval(['p.hair',num2str(h),'Centery = hairdata(3,h+1);']);
    start_id = hairdata(h,1);
    end_id = start_id + hairdata(h,2);
    hair_temp = hair_vertices(start_id:end_id,:);
    eval(['p.h',num2str(h),'_x = hair_temp(:,1);']);
    eval(['p.h',num2str(h),'_y = hair_temp(:,2);']);
    clear hair_temp
end

save([pathbase_data,'hairinfo-files/',Species,...
    '/hairinfo',num2str(run_id),'.mat'],'p')
function hd = modify_hair_data(hd,numofhairs)

for i = 1:numofhairs
    
    hairs(i).x = eval(['hd.hair' num2str(i) 'Centerx']); 
    hairs(i).y = eval(['hd.hair' num2str(i) 'Centery']);
    
end

hd.hairs = hairs; 
hd.radius_m = hd.hdia/2; 
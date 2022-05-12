function [simulation] = concentration_absorbed_by_hairs(simulation)

%global ptindex_hairs hairs_c
%global c 
 
for i=1:length(simulation.inhair_ids)
    %tracking how much concentration is absorbed by each point in each hair
    simulation.hairs_c{i} = simulation.hairs_c{i} + simulation.c(simulation.inhair_ids{i});
    %setting the concentration inside the hairs = 0
    simulation.c(simulation.inhair_ids{i}) = 0; 
end

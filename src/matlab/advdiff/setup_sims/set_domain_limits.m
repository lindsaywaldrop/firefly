function [domainlimits] = set_domain_limits(paths,parameters, run_id)

shift1 = 0.075*parameters.L;
shift2 = 0.02*parameters.L;
runid = str2double(run_id);
hair_vertices = dlmread(strcat(paths.pathbase_data,'vertex-files/',parameters.Species,...
        '/',parameters.Species,'_', num2str(runid),'.vertex'),' ');

far_right = max(max(hair_vertices(2:end,1)))+shift1;
far_left = min(min(hair_vertices(2:end,1)))-shift2;
far_top = max(max(hair_vertices(2:end,2)))+shift2;
far_bottom = 5e-6;

domainlimits = [far_left far_right far_bottom far_top];

%parameters.domainlimits = [-0.5*parameters.L 0.5*parameters.L 0+5e-6 0.35*parameters.L];
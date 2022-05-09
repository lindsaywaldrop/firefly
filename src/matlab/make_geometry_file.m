clear all; close all
%%%% Parameters %%%%
structure_name = "Lucidota_sp";
% Measurements from morphology on antennomere 6
width_olf_hair = 9.67e-6;         % mean width of an olfactory sensillum
width_mech_hair_prox = 5.16e-6;   % mean width of an mechanosensory sensillum (proximal)
width_mech_hair_med = 2.82e-6;    % mean width of an mechanosensory sensillum (medial)
width_mech_hair_dis = 1.28e-6;    % mean width of an mechanosensory sensillum (distal)
olf_hair_length = 13.6e-6;        % Mean length of an olfactory sensillum
mech_hair_length = 77.0e-6;       % Mean length of a mechanosensory sensillum
mech_hair_angle = (29/360)*2*pi;  % Hair angle of a mechanosensory sensillum
mean_dist_all = 25.5e-6;          % Mean distance between any two sensilla in the array
sd_dist_all = 7e-6;               % Standard deviation of distances between sensilla
frac_olf = 0.59;                  % Fraction of all sensilla that are olfactory
density_all = 0.0026/(1e-6)^2;    % Mean density of all sensilla 
width_ant = 490e-6;               % width of antennomere 6

% Parameters for IBM model: 
Nx =  256;       % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Ny =  256;       % # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
Lx = 400e-6;        % Length of Eulerian Grid in x-Direction
Ly = 400e-6;        % Length of Eulerian Grid in y-Direction
dx = Lx/Nx;
dy = Ly/Ny; 
x = 0:dx:Lx-dx;
y = 0:dy:Ly-dy;
[xx, yy] = meshgrid(x,y);
k_Target = 5e6;

% Parameters for 2D model
length_strip = 10e-6; % um, small strip
area_strip = length_strip*width_ant;
total_hairs = ceil(density_all*area_strip);
num_olf_hairs = ceil(frac_olf*total_hairs);
num_mech_hairs = total_hairs-num_olf_hairs;
overlap = floor((mech_hair_length*cos(mech_hair_angle))/mean_dist_all);
base_width = 10e-6;

%%%% create the initial concentration %%%%

xc = 0.9;
[Concentration]= give_Me_Initial_Concentration(Nx,Ny,Lx,Ly,xc);
figure(1)
pcolor(x, y, Concentration), shading flat


%%%% Create the vertex file %%%%%

hairs = zeros(1,total_hairs);
while(sum(hairs)~=num_olf_hairs)
  hairs = 0 + (1-0).*rand(total_hairs,1);
  for j=1:total_hairs
      if(hairs(j)<0.5) 
          hairs(j) = 1;
      else
          hairs(j) = 0;
      end
  end
end

positions = zeros(1,length(hairs));
positions(1) = 50e-6.*rand(1);
for i=2:length(hairs)
  positions(i) = positions(i-1) + (mean_dist_all + sd_dist_all.*randn(1));
end


% Make base plate
x_grid = 0:dx:Lx;
y_grid = 0:dx:base_width;
[pegs_x, pegs_y] = meshgrid(x_grid, y_grid);
pegs(:,1) = reshape(pegs_x,1,[]);
pegs(:,2) = reshape(pegs_y,1,[]);
corners_x = [0 Lx Lx 0];
corners_y = [0,0,base_width,base_width];
[In_base] = inpolygon(xx,yy,corners_x,corners_y);
Concentration(In_base) = 0;

peg_ind_start = length(pegs)+1;

% Make fixed pegs
length_peg(hairs==1) = olf_hair_length; 
length_peg(hairs==0) = mech_hair_length*0.10;
width_peg(hairs==1) = width_olf_hair;
width_peg(hairs==0) = width_mech_hair_prox;
for i=1:2
    for j=1:total_hairs
        if i==1 && hairs(j)==1
            [peg_in, x_peg,y_peg] = make_peg(dx, xx, yy, [positions(j),...
                base_width], length_peg(j), width_peg(j)); 
            pegs = [pegs; [x_peg, y_peg]];
            Concentration(peg_in) = 0;
            clear x_peg y_peg peg_in
        elseif i==2 && hairs(j)==0
            [peg_in, x_peg,y_peg] = make_peg(dx, xx, yy, [positions(j),...
                base_width], length_peg(j), width_peg(j)); 
            pegs = [pegs; [x_peg, y_peg]];
            Concentration(peg_in) = 0;
            clear x_peg y_peg peg_in
        else
            
        end
    end
    if i==1 
        peg_ind_end = length(pegs);
    end
end

save('peg_indices.mat', 'peg_ind_start','peg_ind_end')

% Make flying pegs
for j=1:overlap
  mean_j = ((1/overlap)*1.5*j*mech_hair_length*sin(mech_hair_angle)+base_width);
  sd_j = 0.05*((1/overlap)*j*mech_hair_length*sin(mech_hair_angle)+base_width);
  fly_above = mean_j + sd_j.*randn(num_mech_hairs,1);
  fly_across = (100e-6).*rand(1);
  for i=2:num_mech_hairs
      fly_across = [fly_across, (fly_across(i-1) + (2.5*mean_dist_all + sd_dist_all.*randn(1)))];
  end
  widths(j==1) = width_mech_hair_med;
  widths(j~=1) = width_mech_hair_dis;
  for k=1:num_mech_hairs
      center_point = [fly_across(k) fly_above(k)];
      [peg_in, x_peg, y_peg] = make_flying_peg(dx, xx, yy, center_point, widths);
      pegs = [pegs; [x_peg, y_peg]];
      Concentration(peg_in) = 0;
      clear x_peg y_peg peg_in
  end
end

figure(1)
hold on
plot(pegs(:,1),pegs(:,2),'k.')
xlim([0, Lx])
ylim([0, Ly])
plot(pegs(peg_ind_start:peg_ind_end,1), pegs(peg_ind_start:peg_ind_end,2),'r.')
hold off
print_Lagrangian_Vertices(pegs(:,1), pegs(:,2), structure_name)

print_Lagrangian_Target_Pts(pegs(:,1), k_Target, structure_name);

% Prints .concentration file!

print_Concentration_Info(Nx,Ny,Concentration,structure_name);

CL_save = zeros(length(pegs(peg_ind_start:peg_ind_end,1)),1);
CL_save = CL_save';
dlmwrite('CL_data.csv', CL_save)

clear all

%%%% Other Functions %%%%

function [In, Xin, Yin] = make_peg(dx, xx, yy, center_point, height, width)

  %x_grid = (center_point(1)-0.5*width-1e-6):dx:(center_point(1)+0.5*width+1e-6);
  %y_grid = (center_point(2)-1e-6):dx:(center_point(2)+height+0.5*width+1e-6);
  %[whole_grid_x,whole_grid_y] = meshgrid(x_grid, y_grid);
   whole_grid_x = xx;
   whole_grid_y = yy;
  
  THETA = linspace(0, pi, 25);
  pts_x = 0.5*width*cos(THETA);
  pts_y = 0.5*width*sin(THETA);
  pts_y = pts_y + height + dx + center_point(2);
  pts_x = pts_x + center_point(1);
  pts_x = [pts_x, center_point(1)-0.5*width, 0.5*width+center_point(1)]; 
  pts_y = [pts_y, center_point(2)+dx, center_point(2)+dx];
  
  In = inpolygon(whole_grid_x, whole_grid_y, pts_x, pts_y);
  Xin = whole_grid_x(In);
  Yin = whole_grid_y(In);

end

function [In, Xin, Yin] = make_flying_peg(dx, xx, yy, center_point, width)
  height = width*1.15;
  %x_grid = (center_point(1)-0.5*width-2e-6):dx:(center_point(1)+0.5*width+2e-6);
  %y_grid = (center_point(2)-0.5*height-2e-6):dx:(center_point(2)+0.5*height+2e-6);
  %[whole_grid_x,whole_grid_y] = meshgrid(x_grid, y_grid);
  whole_grid_x = xx;
  whole_grid_y = yy;
  
  THETA = linspace(0, 2*pi, 50);
  pts_x = width*cos(THETA) + center_point(1);
  pts_y = height*sin(THETA) + center_point(2);
  
  In = inpolygon(whole_grid_x, whole_grid_y, pts_x, pts_y);
  Xin = whole_grid_x(In);
  Yin = whole_grid_y(In);
   
end

function print_Lagrangian_Vertices(xLag,yLag,struct_name)

    N = length(xLag);

    vertex_fid = fopen(strcat(struct_name, '.vertex'), 'w');

    fprintf(vertex_fid, '%d\n', N );

    %Loops over all Lagrangian Pts.
    for s = 1:N
        X_v = xLag(s);
        Y_v = yLag(s);
        fprintf(vertex_fid, '%1.16e %1.16e\n', X_v, Y_v);
    end

    fclose(vertex_fid); 
end

function print_Lagrangian_Target_Pts(xLag,k_Target,struct_name)

    N = length(xLag);

    target_fid = fopen(strcat(struct_name,'.target'), 'w');

    fprintf(target_fid, '%d\n', N );

    %Loops over all Lagrangian Pts.
    for s = 1:N
        fprintf(target_fid, '%d %1.16e\n', s, k_Target);
    end

    fclose(target_fid); 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints initial concentration to a file called
% <struct_name>.concentration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

function print_Concentration_Info(Nx,Ny,C,struct_name)

    con_fid = fopen(strcat(struct_name,'.concentration'), 'w');

    for i=1:Ny
        for j=1:Nx
            fprintf(con_fid, '%1.16e ', C(i,j) );
        end
        fprintf(con_fid,'\n');
    end    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives initial concentration gradient inside channel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C] = give_Me_Initial_Concentration(Nx,Ny,Lx,Ly,xc)

x = linspace(0,1,Nx);
y = linspace(0,1,Ny);
[xx,yy]=meshgrid(x,y);
%xc = 0.3; 
C = exp(-7*(2*(xx - xc)./0.4).^2);

end

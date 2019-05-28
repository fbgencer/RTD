%Plot barrier function

function plot_regions(region_matrix)
%size_vec = [x1 y1 width height] 
%rectangle('Position',size_vec,'FaceColor',[1 1 1],'EdgeColor','r','LineWidth',1)
barrier_color = 'k';
wall_thickness = 1;

for iter = 1:size(region_matrix,3)
x1 = region_matrix(:,1,iter);
y1 = region_matrix(:,2,iter);
w = region_matrix(:,3,iter);
h = region_matrix(:,4,iter);

line([x1 x1],[y1 y1+h],'Color',barrier_color,'LineWidth',wall_thickness);
line([x1 x1+w],[y1+h y1+h],'Color',barrier_color,'LineWidth',wall_thickness);
line([x1+w x1+w],[y1 y1+h],'Color',barrier_color,'LineWidth',wall_thickness);
end







res = 0;
end   
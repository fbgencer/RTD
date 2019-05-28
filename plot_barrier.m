%Plot barrier function

function res = plot_barrier(size_vec)
%size_vec = [x1 y1 width height] 
%rectangle('Position',size_vec,'FaceColor',[1 1 1],'EdgeColor','r','LineWidth',1)

x1 = size_vec(1);
y1 = size_vec(2);
w = size_vec(3);
h = size_vec(4);

barrier_color = 'r';
wall_thickness = 1.5;

line([x1 x1],[y1 y1+h],'Color',barrier_color,'LineWidth',wall_thickness);
line([x1 x1+w],[y1+h y1+h],'Color',barrier_color,'LineWidth',wall_thickness);
line([x1+w x1+w],[y1 y1+h],'Color',barrier_color,'LineWidth',wall_thickness);



res = 0;
end   
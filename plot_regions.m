%Plot barrier function

function plot_regions(region_matrix,ax)
%size_vec = [x1 y1 width height] 
%rectangle('Position',size_vec,'FaceColor',[1 1 1],'EdgeColor','r','LineWidth',1)
barrier_color = 'b';
wall_thickness = 2;

nm = 1e-9;
eV = 1.6*10^-19;

x = [];
y = [];

for iter = 1:size(region_matrix,3)
x1 = region_matrix(:,1,iter)/nm;
y1 = region_matrix(:,2,iter)/nm;
w = region_matrix(:,3,iter)/nm;
h = region_matrix(:,4,iter)/eV;
 


%line([x1 x1],[y1 y1+h],'Color',barrier_color,'LineWidth',wall_thickness+1);%sol cizgi
%hold on
xx = linspace(x1,x1+w,10);
yy = (y1+h)*ones(1,size(xx,2));
x = [x,xx];
y = [y,yy];
%xx = linspace(x1,x1+w,10);
%yy = (y1+h)*ones(1,size(x,2));
%plot(xx,yy,'.','Color',barrier_color,'LineWidth',wall_thickness);
%xlabel('nm');
%ylabel('V_o');
%hold on
%line([x1 x1+w],[y1+h y1+h],'Color',barrier_color,'LineWidth',wall_thickness);
%line([x1+w x1+w],[y1 y1+h],'Color',barrier_color,'LineWidth',wall_thickness+1);
end

plot(ax,x,y,'-','Color',barrier_color,'LineWidth',wall_thickness+1);

% interface_x = [0,interface_x,19*nm];
% 
% for iter = 1:size(k,2)
%     xx = linspace(interface_x(iter),interface_x(iter+1),50);
%     psi(:,:,iter) = 0.25+0.01*exp(1i*conj(k(iter))*xx)+0.01*exp(-1i*conj(k(iter))*xx);
%     plot(xx/nm,psi(:,:,iter),'-');
%     hold on    
% end

res = 0;
end   
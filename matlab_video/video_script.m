
clear

load('matlab_video_data_100.mat')

x_hist_t_correct = output_to_matlab.x_hist_t_correct;
new_time_vec = output_to_matlab.new_time_vec;



%%
a_el_1 = [0 90];
a_el_2 = [-90 0];


f1 = figure('units','normalized','outerposition',[0 0 1 1]);

iter = 0;

len = length(1:200:length(x_hist_t_correct));

v = VideoWriter('newfile.avi');
open(v);

[x,y,z] = sphere(40);


% views stuff 
v1 = [-180 30 ];
v2 = [-90 0 ];
v3 = [0 90];
total_T = new_time_vec(end);
vs = [v1' v2' v3'];
ts = [0 total_T/2 total_T];
new_vs = spline(ts,vs,new_time_vec);

for i = 1:10:length(x_hist_t_correct)
% for i = length(x_hist_t_correct)
    cla
    iter = iter+1;
%     iter = length(x_hist_t_correct)
    cla
    figure(f1)
    cla
%     if i ==1 
%         b = axes;
%         b.Color = [.8 .8 .8];
%         cla
%     end
    grid on
    title(['$\Delta$ T = ',num2str(floor(new_time_vec(i))),' Days'],'Interpreter','latex','FontSize',18)
    hold on 
    cla
    set(gcf,'color','w');
%     b = axes;b.Color = 'k';
    hold on
    
    plot3(x_hist_t_correct(1,1:i),x_hist_t_correct(2,1:i),x_hist_t_correct(3,1:i))
%     surf(6378*x,6378*y,6378*z,'FaceColor','blue','EdgeColor', 'black')
    imgRGB = imread('earth.jpg');
    warp(6378*x,6378*y,6378*z,rot90(imgRGB,2))

    plot3(x_hist_t_correct(1,i),x_hist_t_correct(2,i),x_hist_t_correct(3,i),'r.','MarkerSize',28,'LineWidth',4)
    axis equal
    axis(10000*[-6.5 4.5 -4.5 4.5 -1.8 1.8])
%     view((i/50),10)
%     view(-60,5)
    view(new_vs(1,i),new_vs(2,i))
%     title('GTO-GEO Trajectory Optimization')
    xlabel('ECI X (km)','Interpreter','latex','FontSize',18)
    ylabel('ECI Y (km)','Interpreter','latex','FontSize',18)
    zlabel('ECI Z (km)','Interpreter','latex','FontSize',18)
    %view(0 - 90*i/(50*len),90 - 90*i/(50*len))

%     cla
    drawnow
%     set(gcf,'color','w');
%     b = axes;b.Color = 'k';
%     drawnow
     frame= getframe(gcf);
     writeVideo(v,frame)
%      cla(b)
 
    
    
end

close(v);



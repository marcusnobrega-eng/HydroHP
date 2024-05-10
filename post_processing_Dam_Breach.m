%% Post Processing DamBreach-Like
% Developer: Marcus Nobrega
% 9/12/2023


%% Let's Calculate the Forces Associated with the Flow
h_person = 1.70; % m
b_person = 0.35; % m
% gamma_fluid = 1000*9.81; % kg/m3
max_velocity = max(max(Velocity(~isinf(Velocity))));
h = Depth;
area_person = b_person*max(h,h_person);

% Pressure due to velocity
p_vel = (gamma_fluid/g)*Velocity.^2/(2*g); % N/m2

% Force due to velocity
F_vel = p_vel.*area_person/1000/10; % tf

% Pressure due to water depth
% p_depth = (gamma_fluid*h + gamma_fluid*(max(h - h_person,0)))/2*(max(h,h_person)); % N/m2
p_depth = gamma_fluid*h.^2/2 - gamma_fluid*(max(h - h_person,0)).^2/2; % N/m or N per unit width of a person

% Force due to pressuregamma_fluid
F_pressure = p_depth.*b_person/1000/10; % tf
            % N/m * m / 1000 / 10
% Total Force
F_tot = F_vel + F_pressure;
F_tot = F_tot';
max_force = max(max(F_tot));
for i = 1:length(labels.obs_points.nodes)
    node_obs = labels.obs_points.nodes(i);
    max_force_node(i,1) = max(F_tot(node_obs,:)); % tf (force at each observed node)
end

% Water Surface Elevation
wse = Depth + repmat(inv_el',[size(Depth,1),1]);
wse = wse';


%%  Maximum Extents of Forces, Velocities and Depths
color_plot = [21, 179, 196]/255;
color_velocity = [255,99,71]/256;
color_force = [25,25,112]/256;
set(gcf,'units','inches','position',[0,0,7,12])
% Converting x to kilometers
x_plot = x/1000;
max_wse = max(wse')';
label_plot = strcat(labels.simulation_info.ID,'_',labels.simulation_info.NAME,'_','Maximum_WSE_Velocity_Depth.pdf');
    subplot(3,1,1)
    plot(x_plot,inv_el,'LineWidth',4,'LineStyle','-','Color','black')
    hold on
    plot(x_plot,max_wse,'LineWidth',2,'LineStyle','-','Color',color_plot)
    fill([x_plot' fliplr(x_plot')], [inv_el' fliplr(max_wse')],color_plot)
    xlabel('Distance from the dam [km]','Interpreter','latex');
    ylabel('Water surface Elevation [m]','Interpreter','latex');
    ylim([0.98*min(min(wse - Depth')) max(max(1.01*wse))])
    grid on
    title('Maximum Water Surface Elevation','Interpreter','Latex','Fontsize',12);
    hold on
    % Plotting Positions
    for i = 1:length(labels.obs_points.nodes)
        node_obs = labels.obs_points.nodes(i);
        x_node = [node_obs*dx, node_obs*dx]/1000;
        y_node = [inv_el(node_obs), max(max(wse))];
        plot(x_node,y_node,'LineWidth',2,'LineStyle','--','Color','black') 
        hold on
    end
    hold off

    subplot(3,1,2)
    zzz = Velocity;
    zzz(isnan(zzz)) = 0;
    zzz(isinf(zzz)) = 0;
%     velocity(isnan(velocity)) = 0;
    zzz = max(zzz);

    zero_ground = inv_el*0;
    plot(x_plot,zzz,'LineWidth',2,'LineStyle','-','Color',color_velocity)
    xlabel('Distance from the dam [km]','Interpreter','latex');
    ylabel('Velocity [$\mathrm{m~s^{-1}}$]','Interpreter','latex');
    ylim([0 max_velocity])
    grid on
    title('Maximum Velocity','Interpreter','Latex','Fontsize',12);
    hold on
    % Plotting Positions
    for i = 1:length(labels.obs_points.nodes)
        node_obs = labels.obs_points.nodes(i);
        x_node = [node_obs*dx, node_obs*dx]/1000;
        y_node = [0, max(max(Velocity))];
        plot(x_node,y_node,'LineWidth',2,'LineStyle','--','Color','black') 
        hold on
    end
    hold off

    subplot(3,1,3)
    zzz = F_tot';
    zzz(isnan(zzz)) = 0;
    zzz(isinf(zzz)) = 0;
    zzz = max(zzz);
%     zzz(isnan(zzz)) = 0;
    plot(x_plot,zzz,'LineWidth',2,'LineStyle','-','Color',color_force)    
    xlabel('Distance from the dam [km]','Interpreter','latex');
    ylabel('Force [tf]','Interpreter','latex');
    ylim([0 max_force]);
    grid on
    title('Force [tf]','Interpreter','Latex','Fontsize',12);
    hold on
    % Plotting Positions
    for i = 1:length(labels.obs_points.nodes)
        node_obs = labels.obs_points.nodes(i);
        x_node = [node_obs*dx, node_obs*dx]/1000;
        y_node = [0, max(max(F_tot))];
        plot(x_node,y_node,'LineWidth',2,'LineStyle','--','Color','black') 
        hold on
    end
    exportgraphics(gcf,fullfile(folderName,label_plot),'ContentType','vector')   
%%  Video
label_plot = strcat(labels.simulation_info.ID,'_',labels.simulation_info.NAME,'_','Dam_Breach_Like_.mp4');
obj = VideoWriter(fullfile(folderName,label_plot),'MPEG-4');
obj.Quality = 100;
obj.FrameRate = 20;
color_plot = [21, 179, 196]/255;
color_velocity = [255,99,71]/256;
color_force = [25,25,112]/256;
set(gcf,'units','inches','position',[0,0,7,12])
open(obj)
% Converting x to kilometers
x_plot = x/1000;
for n=1:length(time_save)
    hold off
    subplot(3,1,1)
    if n == 1
        t = 1;
        pos = 1;
    else
        t=time_save(n);
        pos = n;
    end
    Plot_Title = 'Time = %.2f (min)';
    sgtitle(sprintf(Plot_Title,t/60),'fontsize',18,'interpreter','latex')    
    plot(x_plot,inv_el,'LineWidth',4,'LineStyle','-','Color','black')
    hold on
    plot(x_plot,wse(:,pos),'LineWidth',2,'LineStyle','-','Color',color_plot)
    fill([x_plot' fliplr(x_plot')], [inv_el' fliplr(wse(:,pos)')],color_plot)
    xlabel('Distance from the dam [km]','Interpreter','latex');
    ylabel('Water surface Elevation [m]','Interpreter','latex');
    ylim([0.98*min(min(wse - Depth')) max(max(1.01*wse))])
    grid on
    title('Water Surface Elevation [m]','Interpreter','Latex');
    hold on
    % Plotting Positions
    for i = 1:length(labels.obs_points.nodes)
        node_obs = labels.obs_points.nodes(i);
        x_node = [node_obs*dx, node_obs*dx]/1000;
        y_node = [inv_el(node_obs), max(max(wse))];
        plot(x_node,y_node,'LineWidth',2,'LineStyle','--','Color','black') 
        hold on
    end
    hold off

    subplot(3,1,2)
    velocity = Velocity(pos,:)';
    velocity(isnan(velocity)) = 0;

    zero_ground = inv_el*0;
    plot(x_plot,velocity,'LineWidth',2,'LineStyle','-','Color',color_velocity)
    fill([x_plot' fliplr(x_plot')], [zero_ground' fliplr(velocity')],color_velocity)
    xlabel('Distance from the dam [km]','Interpreter','latex');
    ylabel('Velocity [m/s]','Interpreter','latex');
    ylim([0 max_velocity])
    grid on
    title('Velocity [m/s]','Interpreter','Latex');
    hold on
    % Plotting Positions
    for i = 1:length(labels.obs_points.nodes)
        node_obs = labels.obs_points.nodes(i);
        x_node = [node_obs*dx, node_obs*dx]/1000;
        y_node = [0, max(max(Velocity))];
        plot(x_node,y_node,'LineWidth',2,'LineStyle','--','Color','black') 
        hold on
    end
    hold off

    subplot(3,1,3)
    zzz = F_tot(:,pos);
    zzz(isnan(zzz)) = 0;
    plot(x_plot,zzz,'LineWidth',2,'LineStyle','-','Color',color_force)
    fill([x_plot' fliplr(x_plot')], [zero_ground' fliplr(zzz')],color_force)
    
    xlabel('Distance from the dam [km]','Interpreter','latex');
    ylabel('Force [tf]','Interpreter','latex');
    ylim([0 max_force]);
    grid on
    title('Force [tf]','Interpreter','Latex');
    hold on
    % Plotting Positions
    for i = 1:length(labels.obs_points.nodes)
        node_obs = labels.obs_points.nodes(i);
        x_node = [node_obs*dx, node_obs*dx]/1000;
        y_node = [0, max(max(F_tot))];
        plot(x_node,y_node,'LineWidth',2,'LineStyle','--','Color','black') 
        hold on
    end

    f = getframe(gcf);
    writeVideo(obj,f);
    hold off    
end
obj.close();


close all

% Video_Name = 'WSE_Top_Width.avi';
% 
% % Set up video
% video = VideoWriter(Video_Name,'MPEG-4');
% open(video);
% 
% % Define water depths for each time
% depths = Depth(:,1)';
% 
% % Preallocate Top Width
% B2 = zeros(size(Flow_Area));
% 
% % Time
% t = time_save; % Sec
% 
% % Define tick size
% ticksize = [0.02 0.01];
% 
% 
% % Color
% color_plot = [21, 179, 196]/255;
% set(gcf,'units','inches','position',[2,0,8,10])
% 
% % Iterate through all time steps
% set(gca,'FontSize',14,'FontName','Garamond')
% for i=1:length(t)
%     Plot_Title = 'Time = %d (sec)';
%     sgtitle(sprintf(Plot_Title, time_store(i)),'fontsize',18,'interpreter','latex')    
%         % --------------- Plotting Channel Width ----------------- %
%     subplot(2,3,[1 2 3]);
%     if flag_section ~= 4
%         B2 = B_function(D,Z1,Z2,a,b,Depth);
%     else
%         for pos_b = 1:length(Flow_Area(1,:))
%             % [y_table, A, P, Rh, y_bar, n_med, Beta, v, B, Q]
%             % [   1,    2, 3, 4,    5,    6,      7,  8,  9, 10]
%             B2(i,pos_b) = Vlookup_g(irr_table,col1,Flow_Area(i,pos_b),9);
%         end
%     end
%     if flag_section == 1
%          offset = b/2 + (Z1 + Z2)/2*max(max(depths));
%          xmax_plot = (Z1 + Z2)*max(max(depths)) + b;
%     elseif flag_section == 2
%         offset = D/2;
%         xmax_plot = D;
%     elseif flag_section == 3
%         offset = xmax/2;
%         xmax_plot = xmax;
%     else
%         offset = max(x_cross)/2; % From station data
%         xmax_plot = max(x_cross);
%     end
%     right_margin = B2(i,:)/2 + offset; left_margin = -B2(i,:)/2 + offset;
%     plot(x_plot,right_margin,'k','LineWidth',2); set(gca,'YDir','reverse');
%     hold on
%     plot(x_plot,left_margin,'k','LineWidth',2); set(gca,'YDir','reverse');
%     hold on
%     fill([x_plot' fliplr(x_plot')], [left_margin fliplr(right_margin)],color_plot)
%     xlabel('$x$ [m]','Interpreter','latex');
%     ylabel('Station [m]','Interpreter','latex');
%     ylim([0, xmax_plot]);
%     xlim([0, max(x)]);
%     grid on
%     title(sprintf('$B_{{max}}(t)$ = %.2f m', max(right_margin - left_margin)),'fontsize',16,'interpreter','latex');
%     set(gca,'FontSize',12,'FontName','Garamond')
%     
%     % ---------------------  Ploting Water Surface Elevation ------- %
%     subplot(2,3,[4 5 6])
%     plot(x_plot,inv_el,'LineWidth',4,'LineStyle','-','Color','k');
%     hold on
%     plot(x_plot,wse(i,:),'k','LineWidth',2,'LineStyle','-','Color',color_plot);
%     fill([x_plot' fliplr(x_plot')], [inv_el' fliplr(wse(i,:))],color_plot)
%     xlabel('$x$ [m]','Interpreter','latex');
%     ylabel('Water Surface Elevation [m]','Interpreter','latex');
%     ylim([0.98*min(min(wse - Depth)) max(max(1.01*wse))])
%     grid on
%     title(sprintf('$WSE_{{max}}(t)$ = %.2f m', max(wse(i,:))),'fontsize',16,'interpreter','latex');
% 
%     % Save the frame for the video
%     set(gca,'FontSize',12,'FontName','Garamond')
%     % Set background color and write to video
%     frame = getframe(gcf);
%     writeVideo(video,frame);
%     hold off
% end
% % Close video writer
% close(video);
% close all





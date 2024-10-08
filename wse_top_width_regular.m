% Post-Processing Routine
% Model: HyPro-SWE
% Developer: Marcus Nobrega
% Last Update: 4/29/2023
% Goal: Create animations of WSE and Top Width for regular sections

close all
label_plot = strcat(labels.simulation_info.ID,'_',labels.simulation_info.NAME,'_','WSE_Top_Width.avi');

% Set up video
video = VideoWriter(label_plot,'MPEG-4');
open(video);

% Define water depths for each time
depths = Depth(:,1)';

% Preallocate Top Width
B2 = zeros(size(Flow_Area));

% Time
t = time_save; % Sec

% Define tick size
ticksize = [0.02 0.01];



% Water Surface Elevation
wse = Depth + repmat(inv_el',[size(Depth,1),1]);

% Color
color_plot = [21, 179, 196]/255;
set(gcf,'units','inches','position',[2,0,8,10])

% Iterate through all time steps
set(gca,'FontSize',14,'FontName','Garamond')
for i=1:length(t)
    Plot_Title = 'Time = %d (sec)';
    sgtitle(sprintf(Plot_Title, time_store(i)),'fontsize',18,'interpreter','latex')    
        % --------------- Plotting Channel Width ----------------- %
    subplot(2,3,[1 2 3]);
    if flag_section ~= 4
        B2 = B_function(D,Z1,Z2,a,b,Depth')';
    else
        for pos_b = 1:length(Flow_Area(1,:))
            % [y_table, A, P, Rh, y_bar, n_med, Beta, v, B, Q]
            % [   1,    2, 3, 4,    5,    6,      7,  8,  9, 10]
            B2(i,pos_b) = Vlookup_g(irr_table,col1,Flow_Area(i,pos_b),9);
        end
    end
    if flag_section == 1
         offset = max(b)/2 + max((Z1 + Z2))/2.*(max(max(Depth)));
         xmax_plot = max((Z1 + Z2)*max(max(Depth)) + b);
    elseif flag_section == 2
        offset = D/2;
        xmax_plot = D;
    elseif flag_section == 3
        offset = xmax/2;
        xmax_plot = xmax;
    else
        offset = max(x_cross)/2; % From station data
        xmax_plot = max(x_cross);
    end
    right_margin = B2(i,:)/2 + offset'; left_margin = -B2(i,:)/2 + offset';
    plot(x,right_margin,'k','LineWidth',2); set(gca,'YDir','reverse');
    hold on
    plot(x,left_margin,'k','LineWidth',2); set(gca,'YDir','reverse');
    hold on
    fill([x' fliplr(x')], [left_margin fliplr(right_margin)],color_plot)
    xlabel('$x$ [m]','Interpreter','latex');
    ylabel('Station [m]','Interpreter','latex');
    ylim([0, xmax_plot]);
    xlim([0, max(x)]);
    grid on
    title(sprintf('$B_{{max}}(t)$ = %.2f m', max(right_margin - left_margin)),'fontsize',16,'interpreter','latex');
    set(gca,'FontSize',12,'FontName','Garamond')
    
    % ---------------------  Ploting Water Surface Elevation ------- %
    subplot(2,3,[4 5 6])
    plot(x,inv_el,'LineWidth',4,'LineStyle','-','Color','k');
    hold on
    plot(x,wse(i,:),'k','LineWidth',2,'LineStyle','-','Color',color_plot);
    fill([x' fliplr(x')], [inv_el' fliplr(wse(i,:))],color_plot)
    xlabel('$x$ [m]','Interpreter','latex');
    ylabel('Water Surface Elevation [m]','Interpreter','latex');
    ylim([0.98*min(min(wse - Depth)) max(max(1.01*wse))])
    grid on
    title(sprintf('$WSE_{{max}}(t)$ = %.2f m', max(wse(i,:))),'fontsize',16,'interpreter','latex');

    % Save the frame for the video
    set(gca,'FontSize',12,'FontName','Garamond')
    % Set background color and write to video
    frame = getframe(gcf);
    writeVideo(video,frame);
    hold off
end
% Close video writer
close(video);
close all

%% Maximum Flood Width in the Channel
        % --------------- Plotting Channel Width ----------------- %
    set(gcf,'units','inches','position',[2,0,6.5,4])        
    linetypes = {'-','--','-.',':','-'};
    linewidths = [2 1.75 1.5 1.25 1];
    color_plot_width = [119,136,153;]/255;
    if flag_section ~= 4
        B2 = B_function(D,Z1,Z2,a,b,Depth')';
    else
        for pos_b = 1:length(Flow_Area(1,:))
            % [y_table, A, P, Rh, y_bar, n_med, Beta, v, B, Q]
            % [   1,    2, 3, 4,    5,    6,      7,  8,  9, 10]
            B2(i,pos_b) = Vlookup_g(irr_table,col1,Flow_Area(i,pos_b),9);
        end
    end
    % Calculating the Flood extent for time-steps defined below:
    time_steps_flood_extent = [tf/5 2*tf/5 3*tf/5 4*tf/5 tf]; % minutes
    dt_flood_extent_store = ceil(time_steps_flood_extent*60/time_store(2));
    dt_flood_extent_store(end) = (dt_flood_extent_store(end))-1;
    for i = 1:length(dt_flood_extent_store)
        flood_extent(:,i) = B2(dt_flood_extent_store(i),:);
        legend_plot(i) = strcat('$ t = $',' ',string(time_steps_flood_extent(i)),' min');
    end
    if flag_section == 1
         offset = max(b)/2 + max((Z1 + Z2))/2.*(max(max(Depth)));
         xmax_plot = max((Z1 + Z2)*max(max(Depth)) + b);
    elseif flag_section == 2
        offset = D/2;
        xmax_plot = D;
    elseif flag_section == 3
        offset = xmax/2;
        xmax_plot = xmax;
    else
        offset = max(x_cross)/2; % From station data
        xmax_plot = max(x_cross);
    end
    for i = 1:length(dt_flood_extent_store)
        right_margin = flood_extent(:,i)/2 + offset'; left_margin = -flood_extent(:,i)/2 + offset';
        ax(i) = plot(x,right_margin,'color',color_plot_width,'LineWidth',linewidths(i),'Linestyle',linetypes(i)); set(gca,'YDir','reverse');
        hold on
        plot(x,left_margin,'color',color_plot_width,'LineWidth',linewidths(i),'Linestyle',linetypes(i)); set(gca,'YDir','reverse');
        hold on
    end
        xlabel('$x$ [m]','Interpreter','latex');
        ylabel('Station [m]','Interpreter','latex');
        ylim([0, xmax_plot]);
        xlim([0, max(x)]);
        grid on
        % Plotting Observed Points
         for i = 1:size(labels.obs_points.length,1)
            node_obs = labels.obs_points.nodes(i);
            x_node = [node_obs*dx, node_obs*dx];
            y_node = [0 xmax_plot];
            plot(x_node,y_node,'LineStyle',line_style(i),'LineWidth',line_width(i),'Color','k')
            hold on;
         end          
        set(gca,'FontSize',12,'FontName','Garamond') 
        set(gca,'TickLength',[0.02 0.01])
        set(gca,'TickDir','out')          
        legend(ax,legend_plot,'interpreter','latex');

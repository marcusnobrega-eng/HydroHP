% Post-Processing Routine
% Model: HyPro-SWE
% Developer: Marcus Nobrega
% Last Update: 4/29/2023
% Goal: Create animations of water depth, top width, and water surface
% elevation

close all
close(video);
label_plot = strcat(labels.simulation_info.ID,'_',labels.simulation_info.NAME,'_','Depth_WSE_Top_Width.mp4');

% Set up video
video = VideoWriter(fullfile(folderName,label_plot),'MPEG-4');
open(video);

% Define water depths for each time
depths = Depth(:,1)';

% Preallocate Top Width
B2 = zeros(size(Flow_Area));

% Time
t = time_save; % Sec

% Define tick size
ticksize = [0.015 0.01];

% Define Tick Position
tickposition = 'in';


% Define polygon for the cross-section
polygon = polyshape(x_cross,y_cross);

% Water Surface Elevation
wse = Depth + repmat(inv_el',[size(Depth,1),1]);

% Color
color_plot = [21, 179, 196]/255;
set(gcf,'units','inches','position',[2,0,8,10])

if flag_elapsed_time ~= 1
    % Time Calculation
    time_duration = time_save/3600/24 + Date_Begin;
end

% Iterate through all time steps
set(gca,'FontSize',14,'FontName','Garamond')
for i=1:(length(t))

    if flag_elapsed_time == 1
        Plot_Title = 'Time = %d (sec)';
        sgtitle(sprintf(Plot_Title, time_store(i)),'fontsize',18,'interpreter','latex')
    else
        sgtitle(string(time_duration(i)),'fontsize',18,'interpreter','latex');
    end
    for j = 1:3 % 3 Cross-sections
        if j == 1
            sec = 1;
        elseif j == 2
            sec = ceil(Nx/2);
        else
            sec = Nx;
        end
        depths = Depth(i,sec)';
        hold on
        subplot(3,3,(j))
        % Set title with time and water depth
        % Define the water depth for this time step
        depth_line = depths*ones(1,length(x_cross));
        plot(x_cross, y_cross, '-k', 'LineWidth', 2,Marker='*'); hold on
        % Find where depth line intersects cross-section polygon
        [x_intersect, y_intersect] = polyxpoly(x_cross,y_cross,x_cross,depth_line);
        if length(x_intersect)  > 1
            % Finding Inside Values
            idx1 = x_cross >= x_intersect(1);
            idx2 = x_cross <= x_intersect(end);
            idx = logical(idx1.*idx2); % Both cases
            x_pol = [x_intersect(1), x_cross(idx)', x_intersect(end)];
            y_pol = [y_intersect(1), y_cross(idx)', y_intersect(2)];
            hold on
            % If the depth line intersects the polygon, plot it
            if ~isempty(x_intersect) && ~isempty(y_intersect)
                depth_plot = depth_line(1)*ones(size(x_pol));
                fill([x_pol fliplr(x_pol)], [y_pol fliplr(depth_plot)],color_plot)
            else
                error('Call developer')
            end
        end
        box on
        if j == 1
            ylabel('Depth [m]','Interpreter','latex')
        end
        xlabel('Station [m]','Interpreter','latex')
        title(sprintf('x = %0.2f m, h = %0.2f m', round((sec-1)*dx,2), depths),'fontsize',16,'interpreter','latex');
        set(gca,'FontSize',12,'FontName','Garamond')
        % Set Tick Postion and Tick Size
        set(gca,'TickLength',ticksize)
        set(gca,'TickDir',tickposition)        
    end

    % --------------- Plotting Channel Width ----------------- %
    subplot(3,3,[4 5 6]);
    if flag_section ~= 4
        B2 = B_function(D,Z1,Z2,a,b,y);
    else
        for pos_b = 1:length(Flow_Area(1,:))
            % [y_table, A, P, Rh, y_bar, n_med, Beta, v, B, Q]
            % [   1,    2, 3, 4,    5,    6,      7,  8,  9, 10]
            B2(i,pos_b) = Vlookup_g(irr_table,col1,Flow_Area(i,pos_b),9);
        end
    end
    offset = max(x_cross)/2; % From station data
    right_margin = B2(i,:)/2 + offset; left_margin = -B2(i,:)/2 + offset;
    plot(x,right_margin,'k','LineWidth',2); set(gca,'YDir','reverse');
    hold on
    plot(x,left_margin,'k','LineWidth',2); set(gca,'YDir','reverse');
    hold on
    fill([x' fliplr(x')], [left_margin fliplr(right_margin)],color_plot)
    xlabel('$x$ [m]','Interpreter','latex');
    ylabel('Station [m]','Interpreter','latex');
    ylim([0, max(x_cross)]);
    grid on
    title(sprintf('$B_{{max}}(t)$ = %0.2f m', max(right_margin - left_margin)),'fontsize',16,'interpreter','latex');
    set(gca,'FontSize',12,'FontName','Garamond')
    % Set Tick Postion and Tick Size
    set(gca,'TickLength',ticksize)
    set(gca,'TickDir',tickposition)

    % ---------------------  Ploting Water Surface Elevation ------- %
    subplot(3,3,[7 8 9])
    plot(x,inv_el,'LineWidth',4,'LineStyle','-','Color','k');
    hold on
    plot(x,wse(i,:),'k','LineWidth',2,'LineStyle','-','Color',color_plot);
    fill([x' fliplr(x')], [inv_el' fliplr(wse(i,:))],color_plot)
    xlabel('$x$ [m]','Interpreter','latex');
    ylabel('Water Surface Elevation [m]','Interpreter','latex');
    ylim([0.98*min(min(wse - Depth)) max(max(1.01*wse))])
    grid on
    title(sprintf('$WSE_{{max}}(t)$ = %0.2f m', max(wse(i,:))),'fontsize',16,'interpreter','latex');
    set(gca,'FontSize',12,'FontName','Garamond')
    % Set Tick Postion and Tick Size
    set(gca,'TickLength',ticksize)
    set(gca,'TickDir',tickposition)

    % Save the frame for the video
    % Set background color and write to video
    frame = getframe(gcf);
    writeVideo(video,frame);
    hold off
    clf
end
% Close video writer
close(video);
close all
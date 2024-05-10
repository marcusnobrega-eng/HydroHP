%%% --------- HyProSWE Model ----------- %%%
% States Post-Processing Routine
% Developer: Marcus Nobrega Gomes Junior
% 5/1/2023
% Goal: Solution of 1-D SVE for given cross-section functions of Area, Perimeter, and
% top Width
% If you have any issues, please contact me at
% marcusnobrega.engcivil@gmail.com

% Generate video showing water level profile over time
close all
label_plot = strcat(labels.simulation_info.ID,'_',labels.simulation_info.NAME,'_','Model_Results.mp4');

% Set up video
video = VideoWriter(fullfile(folderName,label_plot),'MPEG-4');
open(video);

% Set up HEC-RAS colors
hec_ras_colors = [52/255 85/255 132/255; 0 1 1; 0 128/255 1; 0 255/255 0; 1 1 0; 1 128/255 0; 1 0 0; 128/255 0 128/255];
% Set up figure
fig = figure('Units','inches','position', [2,.5,12,7.5]);
set(fig,'DefaultTextInterpreter','latex')
set(gca,'FontSize',16)


t = time_store; % time vector
h = Depth; % water level matrix
q = Discharge; % flow rate matrix
v = Velocity; % velocity matrix
f = Froude; % Froude number matrix
c = Courant; % Courant number matrix
z = x; % distance matrix
x_label = 'Distance (m)';

if flag_elapsed_time ~= 1
    % Time Calculation
    time_duration = time_save/3600/24 + Date_Begin;
end

% Loop over each time step and plot water level profile
for i = 1:(length(t)-1)
    % Plot water level profile
    for j = 1:6
        
        if flag_elapsed_time == 1
            Plot_Title = 'Time = %d (sec)';
            sgtitle(sprintf(Plot_Title, time_store(i)),'fontsize',18,'interpreter','latex')
        else
            sgtitle(string(time_duration(i)),'fontsize',18,'interpreter','latex');
        end
        subplot(2,3,j)
        if j == 1
            variable = Depth;
            y_label = '$h(t)~\mathrm{[m]}$';
            title_label = 'Depth';
        elseif j == 2
            variable = Discharge;
            y_label = '$q_1(t)~\mathrm{[m^3/s]}$';
            title_label = 'Discharge';
        elseif j == 3
            variable = Flow_Area;
            y_label = '$q_2(t)~\mathrm{[m^2]}$';
            title_label = 'Flow Area';
        elseif j == 4
            variable = Froude;
            y_label = '$F_r(t)~\mathrm{[-]}$';
            title_label = 'Froude';
        elseif j == 5
            variable = Courant;
            y_label = '$C_r(t)~\mathrm{[-]}$';
            title_label = 'Courant';
        elseif j == 6
            variable = Velocity;
            y_label = '$v(t)~\mathrm{[m/s]}$';
            title_label = 'Velocity';
        end
        plot(x, variable(i,:), 'k-', 'LineWidth', 2)
        % Fill channel with color
        hold on
        fill([x' fliplr(x')], [zeros(size(x))' fliplr(variable(i,:))],hec_ras_colors(j,:))
        xlabel(x_label,'Interpreter','latex')
        ylabel(y_label,'Interpreter','latex')
        title(sprintf(title_label, time_store(i)))
        set(gca,'FontSize',16)
        ylim([min(variable(:)) max(variable(:))])
        grid on

        % Add data tip for maximum water depth
        [max_variable, max_depth_idx] = max(variable(i,:));
        dcm_obj = datacursormode(fig);
        set(dcm_obj, 'Enable', 'on')
        set(dcm_obj, 'UpdateFcn', {@(obj,evt) sprintf('Distance: %.2fm, Water Level: %.2fm', x(max_depth_idx), max_variable)})
        hold on
        scatter(x(max_depth_idx), max_variable, 50, 'MarkerFaceColor', hec_ras_colors(j,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1.5)
        set(gca,'FontName','Garamond','FontSize',12,'FontWeight','Bold','LineWidth', 1.5);      
        set(gca,'TickLength',[0.02 0.01])
        set(gca,'TickDir','out')
        hold off
    end
    % Set background color and write to video
    frame = getframe(gcf);
    writeVideo(video,frame);
    hold off
end
% Close video writer
close(video);



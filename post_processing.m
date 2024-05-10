%%% --------- HydroHP Model ----------- %%%
% Post-Processing Routine
% Developer: Marcus Nobrega Gomes Junior
% 5/1/2023
% Goal: Solution of 1-D SVE for given cross-section functions of Area, Perimeter, and
% top Width
% If you have any issues, please contact me at
% marcusnobrega.engcivil@gmail.com

%% Creating Modeling Results Folder
% Create the folder name
folderName = 'Modeling_Results';

% Check if the folder already exists
if ~exist(folderName, 'dir')
    % If it doesn't exist, create the folder
    mkdir(folderName);
    disp('Folder "Modeling_Results" created successfully!');
else
    disp('Data sucessfully exported in Modeling_Results Folder');
end

%% Post Processing Graphs
clf
close all


color_plot = [21, 179, 196]/255; % You can change it if you want

% Surfplot
t_save = [0:Nat:tt/dt];
t_save(1,1) = 1;
set(gcf,'units','inches','position',[2,0,8,10])
subplot(3,1,1)
surf(x,tint(t_save)/3600,Froude);
view(0,90);
kk = colorbar ; colormap('jet')
kk.TickDirection = 'out';
shading interp
xlabel('Distance from the dam [m]','Interpreter','latex')
ylabel('Elapsed time [h]','Interpreter','latex')
ylabel(kk,'Froude Number','Interpreter','latex')
zlabel ('Froude Number','Interpreter','Latex');
xlim([0 L]);
ylim([0 tt/60/60]);
set(gca,'FontName','Garamond','FontSize',12,'LineWidth', 1.5);
set(gca,'TickLength',[0.02 0.01])
set(gca,'TickDir','out')
box on

subplot(3,1,2)
surf(x,tint(t_save)/60/60,Depth);
view(0,90);
kk = colorbar ; colormap('jet')
kk.TickDirection = 'out';
shading interp
xlabel('Distance from the dam [m]','Interpreter','latex')
ylabel('Elapsed time [h]','Interpreter','latex')
ylabel(kk,'y (m)','Interpreter','latex')
zlabel ('y (m)','Interpreter','Latex');
xlim([0 L]);
ylim([0 tt/60/60]);
set(gca,'FontName','Garamond','FontSize',12,'LineWidth', 1.5);
set(gca,'TickLength',[0.02 0.01])
set(gca,'TickDir','out')
box on

subplot(3,1,3)
wse = Depth + repmat(inv_el',[size(Depth,1),1]);
surf(x,tint(t_save)/60/60,wse);
view(0,90);
kk = colorbar ; colormap('jet')
kk.TickDirection = 'out';
shading interp
xlabel('Distance from the dam [m]','Interpreter','latex')
ylabel('Elapsed time [h]','Interpreter','latex')
ylabel(kk,'WSE (m)','Interpreter','latex')
zlabel ('WSE (m)','Interpreter','Latex');
xlim([0 L]);
ylim([0 tt/60/60]);
set(gca,'FontName','Garamond','FontSize',12,'LineWidth', 1.5);
set(gca,'TickLength',[0.02 0.01])
set(gca,'TickDir','out')
box on
label_plot = strcat(labels.simulation_info.ID,'_',labels.simulation_info.NAME,'_','Surf_Plots.pdf');
try
    exportgraphics(gcf,fullfile(folderName,label_plot),'ContentType','image','Colorspace','rgb','Resolution',600)
catch ME
    warning('PDF not generated. Please try again later.')
end
clf
close all

if flag_section == 2 % circular
% Video
label_plot = strcat(labels.simulation_info.ID,'_',labels.simulation_info.NAME,'_','Circular_Depth.avi');
obj = VideoWriter(fullfile(folderName,label_plot),'Motion JPEG AVI');
obj.Quality = 100;
obj.FrameRate = 20;
open(obj)
set(gcf,'units','inches','position',[2,2,10,3])
    for n=1:1:(Nt/Nat)
        if n == 1
            t = 1;
            pos = 1;
        else
            t=time_save(n);
            pos = n;
        end
        % Circle Function
        xcir = linspace(0,2*pi,100); % 100 points within 0 and 360 deg
        cir = @(r,ctr) [r*cos(xcir)+ctr(1); r*sin(xcir)+ctr(2)];                     
        c1 = cir(D/2, [D/2; D/2]);

        % Boundary Circle
        % (x - xc)^2 + (y - yc)^2 = D^2/4
        % where xc = D/2 and yc = D/2
        xc = D/2; yc = D/2;
        y01 = Depth(pos,1);
        y02 = Depth(pos,ceil(ceil(Nx/2)));
        y03 = Depth(pos,Nx);
        y0_c = [y01; y02; y03];
        % For a given known y, we have to find two xs, such that
        % x^2 + (-2xc)x + ( (y0 - yc)^2 - xc^2 - D^2/4 )
        % or ax^2 + bx + c, with
        % a = 1; b = -2xc; c = (y0 - yc)^2 - xc^2 - D^2/4
        % x = (- b +- sqrt(b^2 - 4ac)) / (2a)
        a = 1;
        b = -2*xc;
        c = xc^2 + (y0_c - yc).^2 - D^2/4;
        Delta = b^2 - 4*a.*c;
        x1 = (-b + sqrt(Delta))/(2*a);
        x2 = (-b - sqrt(Delta))/(2*a);
        % Now we found the intersection of the circle and a line with know
        % depth        
        subplot(1,3,1)
        title(['t = ',num2str(round(round(t,2),0)),' [sec]'])
        ylim([0 D]);
        xlim([0 D]);
        viscircles([D/2 D/2],D/2,'Color','black');        
%         plot(c1(1,:),c1(2,:),'Color','black');
        hold on
        x_water = linspace(x2(1),x1(1),100);
        y_water = repmat(y01,1,100);
        plot(x_water,y_water,'Color',color_plot,'linewidth',2);
%         fill([c1(1,:) fliplr(c1(1,:))], [y_water fliplr(c2(1,:))],color_plot) 
        ylabel('y(m)','Interpreter','latex')
        xlabel('B(m)','Interpreter','latex')
        legend('Entrance','interpreter','latex')        
        hold off
        grid on
        set(gca,'FontName','Garamond','FontSize',12);
        set(gca,'TickLength',[0.02 0.01])
        set(gca,'TickDir','out');
        box on
        % second section        
        subplot(1,3,2)
        title(['t = ',num2str(round(round(t,2),0)),' [sec]'])
        ylim([0 D]);
        xlim([0 D]);        
        viscircles([D/2 D/2],D/2,'Color','black');
        hold on
        x_water = linspace(x2(2),x1(2),100);
        y_water = repmat(y02,1,100);
        plot(x_water,y_water,'Color',color_plot,'linewidth',2);  
        ylabel('y(m)','Interpreter','latex')
        xlabel('B(m)','Interpreter','latex')     
        legend('x = L/2','interpreter','latex')                
        hold off
        legend('L/2','interpreter','latex')        
        % third section      
        grid on
        set(gca,'FontName','Garamond','FontSize',12);
        set(gca,'TickLength',[0.02 0.01])
        set(gca,'TickDir','out');
        box on        
        subplot(1,3,3)
        title(['t = ',num2str(round(round(t/60),0)),' [sec]'])        
        ylim([0 D]);
        xlim([0 D]);        
        viscircles([D/2 D/2],D/2,'Color','black');
        hold on
        x_water = linspace(x2(3),x1(3),100);
        y_water = repmat(y03,1,100);
        plot(x_water,y_water,'color',color_plot,'linewidth',2);             
        hold off         
        ylabel('y(m)','Interpreter','latex')
        xlabel('B(m)','Interpreter','latex')        
        legend('Exit','interpreter','latex')
        grid on
        set(gca,'FontName','Garamond','FontSize',12);
        set(gca,'TickLength',[0.02 0.01])
        set(gca,'TickDir','out');
        box on        
        % Save frame
        title(['t = ',num2str(round(round(t,2),0)),' [sec]'])
        f = getframe(gcf);
        writeVideo(obj,f);
        hold off 
        clf
    end
obj.close();    
end

if flag_section == 3 % paraboloid
% Video
label_plot = strcat(labels.simulation_info.ID,'_',labels.simulation_info.NAME,'_','Parabolic_Depth.avi');
obj = VideoWriter(fullfile(folderName,label_plot),'Motion JPEG AVI');
obj.Quality = 100;
obj.FrameRate = 20;
open(obj)
set(gcf,'units','inches','position',[2,2,10,3])
    for n=1:1:(Nt/Nat)
        if n == 1
            t = 1;
            pos = 1;
        else
            t=time_save(n);
            pos = n;
        end
        % Save frame
        Plot_Title = 'Time = %d (sec)';        
        sgtitle(sprintf(Plot_Title, time_store(n)),'fontsize',18,'interpreter','latex')
        % Parabolic Function
        % y = a*x^2 => xmax = sqrt((ymax/a))
        ymax = max(max(Depth));
        xmax = sqrt(ymax/a); % x to left and right directions
        xpar = linspace(-xmax,xmax,100); % 100 points within -xmax and xmax deg                     
        ypar = a.*xpar.^2;
        % Now we found bottom of the channel
        % We still need to find xleft and xright for a given y
        y01 = Depth(pos,1);
        y02 = Depth(pos,ceil(Nx/2));
        y03 = Depth(pos,Nx);
        y0_c = [y01; y02; y03];
        xright = sqrt(y0_c/a);
        xleft = - xright;
        subplot(1,3,1)
        title(['t = ',num2str(round(round(t/60,2),0)),' [min]'])
        ylim([0 ymax]);
        xlim([0 ymax]);
        plot(xpar,ypar,'Color','black','LineWidth',2);        
        hold on
        x_water = linspace(xleft(1),xright(1),100);
        y_water = linspace(y01,y01,100);
        plot(x_water,y_water,'color',color_plot,'linewidth',2);
        ylabel('y(m)','Interpreter','latex')
        xlabel('B(m)','Interpreter','latex')
        legend('Entrance','interpreter','latex')      
        grid on
        set(gca,'FontName','Garamond','FontSize',12);
        set(gca,'TickLength',[0.02 0.01])
        set(gca,'TickDir','out');        
        hold off
        % second section        
        subplot(1,3,2)
        title(['t = ',num2str(round(round(t/60,2),0)),' [min]'])
        ylim([0 ymax]);
        xlim([0 ymax]);
        plot(xpar,ypar,'Color','black','LineWidth',2);            
        hold on
        x_water = linspace(xleft(2),xright(2),100);
        y_water = linspace(y02,y02,100);
        plot(x_water,y_water,'color',color_plot,'linewidth',2);
        ylabel('y(m)','Interpreter','latex')
        xlabel('B(m)','Interpreter','latex')
        legend('x = L/2','interpreter','latex')     
        grid on
        set(gca,'FontName','Garamond','FontSize',12);
        set(gca,'TickLength',[0.02 0.01])
        set(gca,'TickDir','out');        
        hold off
        % third section        
        subplot(1,3,3)
        title(['t = ',num2str(round(round(t/60,2),0)),' [min]'])        
        ylim([0 ymax]);
        xlim([0 ymax]);
        plot(xpar,ypar,'Color','black','LineWidth',2);            
        hold on
        x_water = linspace(xleft(3),xright(3),100);
        y_water = linspace(y03,y03,100);
        plot(x_water,y_water,'Color',color_plot,'linewidth',2);
        ylabel('y(m)','Interpreter','latex')
        xlabel('B(m)','Interpreter','latex')
        legend('Outlet','interpreter','latex')         
        grid on
        set(gca,'FontName','Garamond','FontSize',12);
        set(gca,'TickLength',[0.02 0.01])
        set(gca,'TickDir','out');        
        f = getframe(gcf);
        writeVideo(obj,f);
        hold off 
        clf
    end
obj.close();    
end


%% Plots
% Time Scale
if flag_elapsed_time == 1
    close all
    flag_date = 1; % 1 min, 2 hour, 3 day, 4 month
    date_string = {'Elapsed time (min)','Elapsed time (h)','Elapsed time (days)','Elapsed time (months)'};
    if flag_date == 1
        time_scale = 1;
    elseif flag_date == 2
        time_scale = 1/60;
    elseif flag_date == 3
        time_scale = 1/60/24;
    else
        time_scale = 1/60/24/30;
    end
    set(gcf,'units','inches','position',[2,0,8,10])   
    line_style = {'-','--','-.',':','-'};
    line_width = [2 1.75 1.5 1 .75];
    subplot(3,2,1)
    for i = 1:size(labels.obs_points.length,1)
        node_obs = labels.obs_points.nodes(i);
        plot(time_save/60*time_scale,Discharge(:,node_obs),'LineStyle',line_style(i),'LineWidth',line_width(i),'Color','k')
        hold on;
    end
    xlabel(date_string(flag_date),'interpreter','latex');
    ylabel('Flow Discharge (m\textsuperscript{3}/s)','Interpreter','latex');
    legend(labels.obs_points.labels,'Interpreter','Latex','location','best')
%     xlim([1,30]);
    subplot(3,2,2)
    % Velocity
     for i = 1:size(labels.obs_points.length,1)
        node_obs = labels.obs_points.nodes(i);
        plot(time_save/60*time_scale,Velocity(:,node_obs),'LineStyle',line_style(i),'LineWidth',line_width(i),'Color','k')
        hold on;
     end
    %%% Normal Depth Velocity %%%
    xlabel(date_string(flag_date),'interpreter','latex');
    ylabel('Velocity (m/s)','Interpreter','latex');
    legend(labels.obs_points.labels,'Interpreter','Latex','location','best')
%     xlim([1,30]);
    subplot(3,2,3)    
    % Water Depth
     for i = 1:size(labels.obs_points.length,1)
        node_obs = labels.obs_points.nodes(i);
        plot(time_save/60*time_scale,Depth(:,node_obs),'LineStyle',line_style(i),'LineWidth',line_width(i),'Color','k')
        hold on;
     end    
    xlabel(date_string(flag_date),'interpreter','latex');
    ylabel('Water Depths (m)','Interpreter','latex');
    legend(labels.obs_points.labels,'Interpreter','Latex','location','best')
%     xlim([1,30]);
    % Froude Number
    subplot(3,2,4)
     for i = 1:size(labels.obs_points.length,1)
        node_obs = labels.obs_points.nodes(i);
        plot(time_save/60*time_scale,Froude(:,node_obs),'LineStyle',line_style(i),'LineWidth',line_width(i),'Color','k')
        hold on;
     end   
    xlabel(date_string(flag_date),'interpreter','latex');
    ylabel('Froude Number','Interpreter','latex');
    legend(labels.obs_points.labels,'Interpreter','Latex','location','best')
%     xlim([1,30]);
    
    % Courant Number
    subplot(3,2,5)
     for i = 1:size(labels.obs_points.length,1)
        node_obs = labels.obs_points.nodes(i);
        plot(time_save/60*time_scale,Froude(:,node_obs),'LineStyle',line_style(i),'LineWidth',line_width(i),'Color','k')
        hold on;
     end  
    xlabel(date_string(flag_date),'interpreter','latex');
    ylabel('Courant Number','Interpreter','latex');
    legend(labels.obs_points.labels,'Interpreter','Latex','location','best')
%     xlim([1,30]);
    
    % Rating Curve
    % Solving for normal Depth
    ymin = min(min(Depth));
    ymax = max(max(Depth));
    hs = 1; % 1 node
    % hs = ceil(1);
    if flag_section ~= 4
            y_m = [ymin:0.01:ymax]'; % meters
            Qn = 1/nm(hs).*A_function(D,Z1(hs),Z2(hs),a,b(hs),y_m).*Rh_function(D,Z1(hs),Z2(hs),a,b(hs),y_m).^(2/3).*I0(hs)^0.5;
    else
        % [y_table, A, P, Rh, y_bar, n_med, Beta, v, B, Q]
        % [   1,    2, 3, 4,    5,    6,      7,  8,  9, 10]
        col1 = 2; % Col with A
        for jj = 1:length(Flow_Area(:,1))
            Qn(jj,1) = Vlookup_g(irr_table,col1,Flow_Area(jj,hs),10); % Attention here
            y_m(jj,1) = Vlookup_g(irr_table,col1,Flow_Area(jj,hs),1);
            rh_i = Vlookup_g(irr_table,col1,Flow_Area(jj,hs),4);
        end
    end
    subplot(3,2,6)
    tbegin = 30; % (steps), considering initial stabilization of the domain
    plot(Discharge(2:end,hs),Depth(2:end,hs),'LineStyle','--','LineWidth',2,'Color','k')
    hold on
    plot(Discharge(2:end,ceil(Nx/2)),Depth(2:end,ceil(Nx/2)),'LineStyle',':','LineWidth',2,'Color','k')
    hold on
    plot(Qn,y_m,'LineStyle','-','LineWidth',2,'Color','k')
    xlabel('Flow Discharge (m\textsuperscript{3}/s)','Interpreter','latex');
    ylabel('Water Depth (m)','Interpreter','latex');
    ylim([ymin 1.1*max([max(y_m),max(y(ceil(Nx)))])]);
    legend('Q(Inlet)','Q(Nx/2)','$Q_{n}$ (L)','Interpreter','Latex','Location','best')
    hold off
    label_plot = strcat(labels.simulation_info.ID,'_',labels.simulation_info.NAME,'_','Summary_Charts.pdf');
    exportgraphics(gcf,fullfile(folderName,label_plot),'ContentType','vector')
    clf
    close all
else
    close all
    % Time Calculation
    time_duration = time_save/3600/24 + Date_Begin;
    set(gcf,'units','inches','position',[2,0,8,10])
    date_string = {''};
    flag_date = 1;
    subplot(3,2,1)
    % Flows
    plot(time_duration,Discharge(:,1),'LineStyle','--','LineWidth',2,'Color','k')
    hold on
    plot(time_duration,Discharge(:,ceil(Nx/2)),'LineStyle',':','LineWidth',2,'Color','k')
    hold on
    plot(time_duration,Discharge(:,Nx),'LineStyle','-','LineWidth',2,'Color','k')
    hold on
    xlabel(date_string(flag_date),'interpreter','latex');
    ylabel('Flow Discharge (m\textsuperscript{3}/s)','Interpreter','latex');
    legend('Entrance','L/2','Outlet','Interpreter','Latex','location','best')
    % Velocity
    subplot(3,2,2)
    plot(time_duration,Velocity(:,1),'LineStyle','--','LineWidth',2,'Color','k')
    hold on
    plot(time_duration,Velocity(:,ceil(Nx/2)),'LineStyle',':','LineWidth',2,'Color','k')
    hold on
    plot(time_duration,Velocity(:,Nx),'LineStyle','-','LineWidth',2,'Color','k')
    %%% Normal Depth Velocity %%%
    xlabel(date_string(flag_date),'interpreter','latex');
    ylabel('Velocity (m/s)','Interpreter','latex');
    legend('Entrance','L/2','Outlet','Interpreter','Latex','Location','best')
    % Water Depth
    subplot(3,2,3)
    plot(time_duration,Depth(:,1),'LineStyle','--','LineWidth',2,'Color','k')
    hold on
    plot(time_duration,Depth(:,ceil(Nx/2)),'LineStyle',':','LineWidth',2,'Color','k')
    hold on
    plot(time_duration,Depth(:,Nx),'LineStyle','-','LineWidth',2,'Color','k')
    xlabel(date_string(flag_date),'interpreter','latex');
    ylabel('Water Depths (m)','Interpreter','latex');
    legend('Entrance','L/2','Outlet','Interpreter','Latex','Location','best')
    % Froude Number
    subplot(3,2,4)
    plot(time_duration,Froude(:,1),'LineStyle','--','LineWidth',2,'Color','k')
    hold on
    plot(time_duration,Froude(:,ceil(Nx/2)),'LineStyle',':','LineWidth',2,'Color','k')
    hold on
    plot(time_duration,Froude(:,Nx),'LineStyle','-','LineWidth',2,'Color','k')
    xlabel(date_string(flag_date),'interpreter','latex');
    ylabel('Froude Number','Interpreter','latex');
    legend('Entrance','L/2','Outlet','Interpreter','Latex','Location','best')
    
    % Courant Number
    subplot(3,2,5)
    plot(time_duration,Courant(:,1),'LineStyle','--','LineWidth',2,'Color','k')
    hold on
    plot(time_duration,Courant(:,ceil(Nx/2)),'LineStyle',':','LineWidth',2,'Color','k')
    hold on
    plot(time_duration,Courant(:,Nx),'LineStyle','-','LineWidth',2,'Color','k')
    xlabel(date_string(flag_date),'interpreter','latex');
    ylabel('Courant Number','Interpreter','latex');
    legend('Entrance','L/2','Outlet','Interpreter','Latex','Location','best')
    
    % Rating Curve
    % Solving for normal Depth
    ymin = min(min(Depth));
    ymax = max(max(Depth));
    hs = 1; % 1 node
    % hs = ceil(1);
    if flag_section ~= 4
            y_m = [ymin:0.01:ymax]'; % meters
            Qn = 1/nm(hs).*A_function(D,Z1,Z2,a,b,y_m).*Rh_function(D,Z1,Z2,a,b,y_m).^(2/3).*I0(hs)^0.5;
    else
        % [y_table, A, P, Rh, y_bar, n_med, Beta, v, B, Q]
        % [   1,    2, 3, 4,    5,    6,      7,  8,  9, 10]
        col1 = 2; % Col with A
        for jj = 1:length(Flow_Area(:,1))
            Qn(jj,1) = Vlookup_g(irr_table,col1,Flow_Area(jj,hs),10); % Attention here
            y_m(jj,1) = Vlookup_g(irr_table,col1,Flow_Area(jj,hs),1);
            rh_i = Vlookup_g(irr_table,col1,Flow_Area(jj,hs),4);
        end
    end
    subplot(3,2,6)
    tbegin = 30; % (steps), considering initial stabilization of the domain
    plot(Discharge(2:end,hs),Depth(2:end,hs),'LineStyle','--','LineWidth',2,'Color','k')
    hold on
    plot(Discharge(2:end,ceil(Nx/2)),Depth(2:end,ceil(Nx/2)),'LineStyle',':','LineWidth',2,'Color','k')
    hold on
    plot(Qn,y_m,'LineStyle','-','LineWidth',2,'Color','k')
    xlabel('Flow Discharge (m\textsuperscript{3}/s)','Interpreter','latex');
    ylabel('Water Depth (m)','Interpreter','latex');
    ylim([ymin 1.1*max([max(y_m),max(y(ceil(Nx)))])]);
    legend('Q(Inlet)','Q(Nx/2)','$Q_{n}$ (L)','Interpreter','Latex','Location','best')
    hold off
    label_plot = strcat(labels.simulation_info.ID,'_',labels.simulation_info.NAME,'_','Summary_Charts.pdf');    
    exportgraphics(gcf,fullfile(folderName,label_plot),'ContentType','vector')
    clf
    close all
end
%% Dam Breach - Like Post Processing
% post_processing_Dam_Breach

%% States Post-Processing
% states_post_processing
%% Cross-Section Post-Processing
% if flag_section == 4
%     cross_section_post_processing
% end
%% Lateral Profiles
% if flag_section ~= 4
%     wse_top_width_regular
% end

%% Detailed Output
Detailed_Output_Script

%% Dam Breach-like Detailed Output
Detailed_Output_Dam_Breach_Like
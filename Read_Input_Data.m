%%% --------- HyProSWE Model ----------- %%%
% Script to read input data
% Developer: Marcus Nobrega Gomes Junior
% 5/1/2023
% Goal: Solution of 1-D SVE for given cross-section functions of Area, Perimeter, and
% top Width
% If you have any issues, please contact me at
% marcusnobrega.engcivil@gmail.com

% ---------- Please, don't change anything below ------------ %


%% Read Input Data %%
data = readtable(HydroHP_Input_File,'Sheet','Input_Data');

data_labels = readtable(HydroHP_Input_File,'Sheet','Input_Data');

% Labels
labels.simulation_info.ID = string((table2array(data_labels(2,17))));
labels.simulation_info.NAME = string(table2cell(data_labels(2,18)));

% Observation Points
labels.obs_points.length = (table2array(data_labels(2:6,22)));
n_points = sum(~isnan(labels.obs_points.length));

labels.obs_points.length = (table2array(data_labels(2:(2 + n_points -1),22)));

labels.obs_points.nodes = (table2array(data_labels(2:(2+n_points-1),20)));

labels_name = ((data_labels(2:(2 + n_points-1),21)));


for i = 1:n_points
    labels.obs_points.labels{i,:} = string(labels_name{i,:});
end


data_2 = xlsread(HydroHP_Input_File,'Input_Data');

b = 0; Z1 = 0; Z2 = 0; a = 0; D = 0;


% General Data
general_data = table2array(data(1:17,2));
L = general_data(1,1);
Nx = general_data(2,1);
el = general_data(3,1);
g = general_data(4,1);
nm = general_data(5,1);
I0 = general_data(6,1);
tf = general_data(7,1);
dt = general_data(8,1);
animation_time = general_data(9,1);
s_outlet = general_data(10,1);
dh = general_data(11,1);
alpha = general_data(12,1);
dtmin = general_data(13,1);
dtmax = general_data(14,1);
gamma_fluid = general_data(15,1);


% Flags
flags = table2array(data(19:32,2));
flag_hydrograph = flags(1,1);
flag_outlet = flags(2,1);
flag_friction = flags(3,1);
flag_section = flags(4,1);
flag_stage_hydrograph = flags(5,1);
flag_nash = flags(6,1);
flag_slope = flags(7,1);
flag_elevation = flags(8,1);
flag_output = flags(9,1);
flag_plot_HP = flags(10,1);
flag_elapsed_time = flags(11,1);
flag_manning = flags(12,1);
flag_trapezoid_sections = flags(13,1);
flag_dam_break_hydrograph = flags(14,1);

% Dam Info
if flag_dam_break_hydrograph == 1
    dam_data = (table2array(data(2:4,25)));  
    h_dam = dam_data(1);
    W_dam = dam_data(2);
    B_dam = dam_data(3);
end


if flag_elapsed_time ~= 1
    Date_Begin = general_data(16,1);
    Date_Begin = datetime(datestr(Date_Begin+datenum('30-Dec-1899')));
    Date_End = general_data(17,1);
    Date_End = datetime(datestr(Date_End+datenum('30-Dec-1899')));    
end

if flag_nash == 1
    nash_data = table2array(data(1:4,5));
    % Hydrograph
    Tp = nash_data(1,1);
    Qb = nash_data(2,1);
    Beta = nash_data(3,1);
    Qp = nash_data(4,1);
else
    % Input Hydrograph
    input_hydrograph_data = table2array(data(8:end,4:5)); % m3/s
    time_ = input_hydrograph_data(1:end,1); % min
    Qe1_ = input_hydrograph_data(1:end,2);
%     [time_,Qe1_] = breach_hydrograph(h,W,B,dt,tf,flag_plot);   
    Qe1 = zeros(size(Qe1_,1) - sum(isnan(Qe1_)),1);
    time = zeros(size(time_,1) - sum(isnan(time_)),1);    
    % Taking away nans
    for i = 1:length(Qe1)
        if isnan(Qe1_(i)) || isnan(time_(i))
            break
        else
            Qe1(i,1) = Qe1_(i,1);
            time(i,1) = time_(i,1);
        end
    end
    clear Qe1_ time_
end

if flag_stage_hydrograph ~= 0
    % Stage Hydrograph
    input_stage_data = table2array(data(8:end,7:8));
    time_stage_ = input_stage_data(1:end,1);
    he1_ = data(1:end,2);
    he1 = zeros(size(he1_,1) - sum(isnan(he1_)),1);
    time_stage = zeros(size(time_stage_,1) - sum(isnan(time_stage_)),1);
    % Taking away nans
    for i = 1:length(he1)
        if isnan(he1_(i)) || isnan(time_stage_(i))
            break
        else
            he1(i,1) = he1_(i,1);
            time_stage(i,1) = time_stage_(i,1);
        end
    end
    clear Qe1_ time_stage_
end

if flag_slope ~= 0
    % Slope
    input_slope_data = table2array(data(8:end,11:12));
    station = input_slope_data(1:end,1);
    bottom_slope = input_slope_data(2:end,2);
    slopes_not_nan = zeros(size(bottom_slope,1) - sum(isnan(bottom_slope)),1);
    station_index = zeros(size(station,1) - sum(isnan(station)),1);
    % Taking away nans
    for i = 1:length(station_index)
        if isnan(bottom_slope(i)) || isnan(station(i))
            break
        else
            slopes_not_nan(i,1) = bottom_slope(i,1);
            station_index(i,1) = station(i,1);
        end
    end
    clear station_index station bottom_slope
    bottom_slope = slopes_not_nan;
end

if flag_elevation ~= 0
    % Slope
    input_slope_data = table2array(data(8:end,7:8));
    station = input_slope_data(1:end,1);
    elevation_cell = table2array(data(8:end,7:8));
    inv_el_ = zeros(size(elevation_cell,1) - sum(isnan(elevation_cell)),1);
    station_index = zeros(size(station,1) - sum(isnan(station)),1);
    % Taking away nans
    for i = 1:length(station_index)
        if isnan(elevation_cell(i)) || isnan(station(i))
            break
        else
            inv_el_(i,1) = elevation_cell(i,1);
            station_index(i,1) = station(i,1);
        end
    end
    inv_el = inv_el_; % Invert Elevation
    clear station_index station elevation_cell inv_el_
end

if flag_manning == 1
    % Manning
    input_manning_data = table2array(data(8:end,14));
    manning_values = zeros(size(input_manning_data,1) - sum(isnan(input_manning_data)),1);
    % Taking away nans
    for i = 1:length(manning_values)
        if isnan(manning_values(i)) || isnan(input_manning_data(i))
            break
        else
            manning_values(i,1) = input_manning_data(i,1);
        end
    end
    clear input_manning_data
end

% Trapezoid Cross-Section Data
if flag_trapezoid_sections == 1

    input_trapezoid_data = table2array(data(8:end,15:17));
    b_values = zeros(size(input_trapezoid_data,1) - sum(isnan(input_trapezoid_data(:,1))),1);
    slope_left_values = zeros(size(input_trapezoid_data,1) - sum(isnan(input_trapezoid_data(:,1))),1);
    slope_right_values = zeros(size(input_trapezoid_data,1) - sum(isnan(input_trapezoid_data(:,1))),1);
    % Taking away nans
    for i = 1:length(b_values)
        if isnan(b_values(i)) || isnan(input_trapezoid_data(i,1))
            break
        else
            b_values(i,1) = input_trapezoid_data(i,1);
            slope_left_values(i,1) = input_trapezoid_data(i,2);
            slope_right_values(i,1) = input_trapezoid_data(i,3);
        end
    end
    % Read all cross-section data
    cross_section.trapezoid.b = b_values;
    cross_section.trapezoid.slope_left = slope_left_values;
    cross_section.trapezoid.slope_right = slope_right_values;
    clear input_trapezoid_data slope_right_values slope_left_values b_values
end

% Outlet
if flag_outlet ~=1
    input_slope_wave = table2array(data(8:5,8));
    h_0_wave = input_slope_wave(1,1);
    H_0_wave = input_slope_wave(2,1);
    L_wave = input_slope_wave(3,1);
    T_wave = input_slope_wave(4,1);
    x_wave = input_slope_wave(5,1);
end

% Section
if flag_section == 1
    input_slope_trapezoid = table2array(data(1:3,11));
    b = (input_slope_trapezoid(1,1))*ones(Nx,1);
    Z1 = input_slope_trapezoid(2,1)*ones(Nx,1);
    Z2 = input_slope_trapezoid(3,1)*ones(Nx,1);
elseif flag_section == 2
    input_slope_circular = table2array(data(1,14));
    D = input_slope_circular(1,1);
elseif flag_section == 3
    input_slope_parabolic = table2array(data(3,14));
    a = data(1,1);
else
    % Read HP estimator data
    [y_irr, A_irr, P_irr, Rh_irr, y_bar_irr, n_med_irr, Beta_irr, u_irr, B_irr, Q_irr, x_cross, y_cross,s0] = HP_estimator(flag_plot_HP,dh);
    irr_table = [y_irr, A_irr, P_irr, Rh_irr, y_bar_irr, n_med_irr, Beta_irr, u_irr, B_irr, Q_irr];

    % Some Boundary Conditions
    % [y_irr, A_irr, P_irr, Rh_irr, y_bar_irr, n_med_irr, Beta_irr, u_irr, B_irr, Q_irr];
    % [   1,    2,     3,      4,    5,            6,          7,  8,      9,  10]
    irr_table(1,6) = irr_table(2,6); irr_table(1,7) = 0; irr_table(1,8) = 0;
    % Second Line
    irr_table(2,2) = 0; irr_table(2,3) = 0; irr_table(2,4) = 0; irr_table(2,5) = 0; irr_table(2,7) = 0; irr_table(2,8) = 0; irr_table(2,9) = 0; 
%     z = irr_table;
%     second = 0*z(2,:);
%     second(1,1) = 0.5*10^-3; second(1,6) = z(1,6); second(1,10) = z(1,10)/2;
%     z = [z(1,:) ; second; z(2:end,:)];
%     irr_table = z;
end


% Contraint at observed flow
if flag_hydrograph == 1
    if max(time) ~= tf
        z = round(tf - max(time),0);
        for i = 1:z
            Qe1(end + 1,1) = 0;
            time(end+1,1) = time(end,1) + 1;
        end
    end
end


% Contraint at stage hydrograph
if flag_stage_hydrograph == 1
    if max(time_stage) ~= tf
        z = round(tf - max(time_stage),0);
        for i = 1:z
            he1(end + 1,1) = 0;
            time_stage(end+1,1) = time_stage(end,1) + 1;
        end
    end
end

%% Normal Input Data

% Manning
dx = L/(Nx-1);
if flag_manning ~= 1
    for i = 1:Nx
        length_downstream = (i-1)*dx;    
        if length_downstream <= Lf + dx
            nm(i,1) = nf;
        else
            nm(i,1) = nc;
        end
    end
else
    % Manning
      nm = manning_values;
end

% Hydrograph
if  flag_dam_break_hydrograph == 1
    flag_plot = 1;
    [time_,Qe1_] = breach_hydrograph(h_dam,W_dam,B_dam,0.1,tf*60,flag_plot,labels.simulation_info.ID,labels.simulation_info.NAME);
    Qe1 = zeros(size(Qe1_,1) - sum(isnan(Qe1_)),1);
    time = zeros(size(time_,1) - sum(isnan(time_)),1);
end

% Boundary condition at t = 0
if flag_dam_break_hydrograph == 1
    Qe1_(1,1) = 5; % 5m3/s as an initial flow hydrograph
    warning('Assuming a 5m3/s initial flow.')
    pause(1.5);
end

% Taking away nans
if flag_dam_break_hydrograph == 1
    for i = 1:length(Qe1)
        if isnan(Qe1_(i)) || isnan(time_(i))
            break
        else
            Qe1(i,1) = Qe1_(i,1);
            time(i,1) = time_(i,1);
        end
    end
end

%% Entering Automatic Input Data

% % Dam IDs
% Dam_IDs = [5	20	21	22	26	30	31	32	33	40	41	54	75	149	432	1097	1910	7112	7120	19610];
% 
% % Dam Names
% Dam_Names = {'Fazenda São Pedro','Caldeirões', 'Cocorobó','Rio Paranã', 'Cacimba da Várzea', 'Escondido I', 'Felismina Queiroz', 'Mãe d Água', 'Poleiros', 'Várzea Grande', 'Carnaúba', 'Estreito', 'Ipanema I', 'Sítio Carnaúba', 'Zabumbão', 'Pedro Beicht', 'Itapajé', 'Espinheiro', 'Todos os Santos', 'Fazenda Manaus'};
% 
% % Dam Information
% Input_Data_Dams = [9.00	178.00	1492.84	8827.385	1743.787	0.035	0.025	100.00	15.06	13.95	0.001825601
% 29.00	297.00	218.00	1000.00	4320.00	0.035	0.025	80.00	2.86	2.15	0.00566
% 33.50	1320.00	5549.00	72.00	257.00	0.035	0.025	1000.00	26.03	26.03	0.00167
% 33.00	1760.00	2926.00	6475.00	550.00	0.035	0.025	5.00	6.39	3.63	0.00731
% 22.30	311.00	1336.00	1000.00	500.00	0.035	0.025	180.00	8.78	8.78	0.00157
% 12.50	1200.00	1109.00	100.00	1470.00	0.035	0.025	340.00	33.69	33.69	0.00096
% 13.00	391.00	330.00	140.00	655.00	0.035	0.025	30.00	12.71	13.30	0.00536
% 35.00	175.00	7837.00	816.00	601.00	0.035	0.025	75.00	15.89	15.89	0.0051
% 25.00	430.00	2653.00	3055.00	3169.00	0.035	0.025	100.00	11.91	26.03	0.00566
% 25.00	544.00	2140.00	5145.25	1633.80	0.035	0.025	150.00	7.60	6.77	0.00196
% 19.00	550.00	3014.00	14600.00	366.00	0.035	0.025	50.00	31.82	30.14	0.00218
% 28.00	1091.00	2483.00	2165.00	1982.00	0.035	0.025	300.00	4.17	30.14	0.00132
% 16.40	724.00	218.00	132.00	404.00	0.035	0.025	10.00	19.08	8.51	0.01829
% 6.00	329.00	250.00	14500.00	1500.00	0.035	0.025	100.00	5.73	8.03	0.04855
% 65.00	365.00	2565.00	1171.00	4974.00	0.035	0.025	250.00	5.34	4.26	0.00167
% 23.00	347.00	2151.00	3711.00	10590.00	0.035	0.025	100.00	9.68	13.00	0.00284
% 17.90	436.00	2193.00	2170.00	685.00	0.035	0.025	1300.00	4.01	5.14	0.00584
% 4.00	115.00	652.00	9010.00	1780.00	0.035	0.025	30.00	9.51	9.51	0.00209
% 31.00	35.00	10808.00	4990.00	5710.00	0.035	0.025	50.00	2.66	2.79	0.00501
% 7.00	156.00	276.00	3896.00	1317.00	0.035	0.025	70.00	38.19	22.02	0.00244];
% 
% % h , B , W, Lf, Lc, nf, nc, b, slope_left, slope_right, I0
% 
% 
% h = Input_Data_Dams(ii,1); % Height of the dam
% B = Input_Data_Dams(ii,2); % Width of the dam
% W = Input_Data_Dams(ii,3); % Length of the dam
% Lf = Input_Data_Dams(ii,4); % Lentgh of the floodplain [m]
% Lc = Input_Data_Dams(ii,5); % Length of the city
% nf = Input_Data_Dams(ii,6); % Manning's roughness of the floodplain
% nc = Input_Data_Dams(ii,7); % Manning's roughness of the city
% b = Input_Data_Dams(ii,8)*ones(Nx,1); % Base of the channel [m]
% slope_left = Input_Data_Dams(ii,9)*ones(Nx,1); % Left cotang of the channel
% slope_right = Input_Data_Dams(ii,10)*ones(Nx,1); % Right cotang of the channel
% I0 = Input_Data_Dams(ii,11); % Average slope
% s_outlet = I0;
% 
% L = Lf + Lc + 1000; % Total Length [m]
% 
% % Manning
% dx = L/(Nx-1);
% for i = 1:Nx
%     length_downstream = (i-1)*dx;    
%     if length_downstream <= Lf + dx
%         nm(i,1) = nf;
%     else
%         nm(i,1) = nc;
%     end
% end
% 
% if flag_manning == 1
%     % Manning
%       nm = manning_values;
% end
% 
% % Labels
% labels.simulation_info.ID = string(Dam_IDs(ii));
% labels.simulation_info.NAME = string(Dam_Names(ii));
% 
% % Hydrograph
% if flag_dam_break_hydrograph == 1
%     flag_plot = 1;
%     [time_,Qe1_] = breach_hydrograph(h,W,B,0.1,tf*60,flag_plot,labels.simulation_info.ID,labels.simulation_info.NAME);
%     Qe1 = zeros(size(Qe1_,1) - sum(isnan(Qe1_)),1);
%     time = zeros(size(time_,1) - sum(isnan(time_)),1);
% end
% 
% % Boundary condition at t = 0
% Qe1_(1,1) = 5; % 5m3/s as an initial flow hydrograph 
% 
% % Taking away nans
% for i = 1:length(Qe1)
%     if isnan(Qe1_(i)) || isnan(time_(i))
%         break
%     else
%         Qe1(i,1) = Qe1_(i,1);
%         time(i,1) = time_(i,1);
%     end
% end
% 
% % Observed Nodes
% labels.obs_points.length = [dx (dx + Lf) (dx + Lf + Lc) (dx + Lf + Lc + 1000)]';
% 
% labels.obs_points.nodes = round((labels.obs_points.length)/dx,0)';




% Dam-Breach Like Outputs

% Outputs
depth_flood = Depth(1,:) + 0.5; % Initial Depth + 5 cm
clear time_max_depth max_depth_node time_max_velocity max_velocity_node time_arrival_node max_force_node computational_time

% Force
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

for i = 1:length(labels.obs_points.nodes)
    node_obs = labels.obs_points.nodes(i);
    position = find(Depth(:,node_obs) == max(Depth(:,node_obs)),1,'first');
    time_max_depth(i) = time_store(position)/60; % min
    max_depth_node(i) = Depth(position,node_obs); % m
    position = find(Velocity(:,node_obs) == max(Velocity(:,node_obs)),1,'first');
    time_max_velocity(i) = time_store(position)/60; % min
    max_velocity_node(i) = Velocity(position,node_obs); % m/s
    position = find(Depth(:,node_obs) > depth_flood(node_obs),1,'first');
    if isempty(position)
        time_arrival_node(i) = nan;
    else
        time_arrival_node(i) = time_store(position)/60; % min
    end
    max_force_node(i,1) = max(F_tot(node_obs,:));
end

% Adding extra values to be 5 entries
n_extra = 5 - length((labels.obs_points.nodes));
size_data = size(time_max_depth,2);
time_max_depth((size_data+1):(size_data+ n_extra)) = nan;
max_depth_node((size_data+1):(size_data+ n_extra)) = nan;
time_max_velocity((size_data+1):(size_data+ n_extra)) = nan;
max_velocity_node((size_data+1):(size_data+ n_extra)) = nan;
time_arrival_node((size_data+1):(size_data+ n_extra)) = nan;
max_force_node((size_data+1):(size_data+ n_extra),1)= nan;
computational_time = toc;
%% Outputs
final_parameters = [computational_time/60, tf/60, dx, ...
    time_max_depth, max_depth_node, time_max_velocity, max_velocity_node, time_arrival_node, ...
    max_force_node'];

%% Table
headers = [{'Computational Time (min)'},{'Simulation Time (min)'},{'Spatial Discretization (m)'},...
    {'Arrival time of maximum depth (min) [1]'}, {'Arrival time of maximum depth (min) [2]'}, {'Arrival time of maximum depth (min) [3]'}, {'Arrival time of maximum depth (min) [4]'}, {'Arrival time of maximum depth (min) [5]'}, ...
    {'Maximum Depth (m) [1]'}, {'Maximum Depth (m) [2]'}, {'Maximum Depth (m) [3]'}, {'Maximum Depth (m) [4]'}, {'Maximum Depth (m) [5]'}, ...
    {'Arrival time of maximum velocity [1]'},{'Arrival time of maximum velocity [2]'}, {'Arrival time of maximum velocity [3]'}, {'Arrival time of maximum velocity [4]'}, {'Arrival time of maximum velocity [5]'}, ...
    {'Maximum velocity (m/s) [1]'}, {'Maximum velocity (m/s) [2]'}, {'Maximum velocity (m/s) [3]'} , {'Maximum velocity (m/s) [4]'}, {'Maximum velocity (m/s) [5]'}, ...
    {'Flood arrival time (min) [1]'}, {'Flood arrival time (min) [2]'}, {'Flood arrival time (min) [3]'}, {'Flood arrival time (min) [4]'}, {'Flood arrival time (min) [5]'}, ...
    {'Maximum force (tf) [1]'}, {'Maximum force (tf) [2]'}, {'Maximum force (tf) [3]'}, {'Maximum force (tf) [4]'}, {'Maximum force (tf) [5]'}];
data_save = [final_parameters]; % Concatenating dataset to the time
T = array2table(data_save);
T.Properties.VariableNames(1:size(data_save,2)) = headers;
label_plot = strcat(labels.simulation_info.ID,'_',labels.simulation_info.NAME,'_','Dam_Breach_Like_Outputs.csv');
writetable(T,fullfile(folderName,label_plot),'Delimiter',',');

disp('Attention: Data exported in .CSV');

zzz = [time_arrival_node(2) max_depth_node(2)  ; time_arrival_node(3) max_depth_node(3); time_arrival_node(4) max_depth_node(4)  ];
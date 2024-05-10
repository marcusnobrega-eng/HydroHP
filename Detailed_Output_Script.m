% HydroHP Model
% Output .csv script
% Developer: Marcus Nobrega
% Goal: Create a detailed output from modeling results
% Last updated: 4/30/2023


%%% ----------------------- All rights reserved --------------------- %%

% Number of states
ns = 6;
% 0 - time, 1 - flow, 2 - depth, 3 - velocity, 4 - Courant, 5 - Froude, 6,
% 7 WSE

% Concatenate data
t = time_store; % time vector
h = Depth; % water level matrix
q = Discharge; % flow rate matrix
v = Velocity; % velocity matrix
f = Froude; % Froude number matrix
c = Courant; % Courant number matrix
z = x; % distance matrix

% Round Data
decimal_places = 3;

% Returning wse to correct size
% wse = wse';

data = zeros(size(Depth,1),size(Depth,2),ns);

data(:,:,1) = Depth;
data(:,:,2) = Discharge;
data(:,:,3) = Velocity;
data(:,:,4) = Froude;
data(:,:,5) = Courant;
data(:,:,6) = wse;


if flag_output == 1
    for i = 1:(Nx*ns)
        j = floor((i-1)/ns);
        x_cell = j*dx;
        if mod(i-1,ns) == 0 || (i-1)/ns == 1
            states_title(1,i) = cellstr(sprintf('Depth (m), x(m) = %0.2f',x_cell));
        elseif mod(i-1,ns) == 1 || (i-1)/ns == 2
            states_title(1,i) = cellstr(sprintf('Discharge (m^3/s), x(m) = %0.2f',x_cell));
        elseif mod(i-1,ns) == 2 || (i-1)/ns == 3
            states_title(1,i) = cellstr(sprintf('Velocity (m/s), x(m) = %0.2f',x_cell));
        elseif mod(i-1,ns) == 3 || (i-1)/ns == 4
            states_title(1,i) = cellstr(sprintf('Froude (-), x(m) = %0.2f',x_cell));
        elseif mod(i-1,ns) == 4 || (i-1)/ns == 5
            states_title(1,i) = cellstr(sprintf('Courant Number (-), x(m) = %0.2f',x_cell));
        elseif mod(i-1,ns) == 5 || (i-1)/ns == 6
            states_title(1,i) = cellstr(sprintf('Water Surface Elevation (m), x(m) = %0.2f',x_cell));
        end
    end
else
    for i = 1:(Nx*ns)
        if mod(i,Nx) ~= 0
            j = mod(i,Nx);
            x_cell = (j-1)*dx; % m
        else
            j = Nx;
            x_cell = (j-1)*dx; % m
        end
        if floor(i/Nx) == 0 || i/Nx == 1
            states_title(1,i) = cellstr(sprintf('Depth (m), x(m) = %0.2f',x_cell));
        elseif floor(i/Nx) == 1 || i/Nx == 2
            states_title(1,i) = cellstr(sprintf('Discharge (m^3/s), x(m) = %0.2f',x_cell));
        elseif floor(i/Nx) == 2 || i/Nx == 3
            states_title(1,i) = cellstr(sprintf('Velocity (m/s), x(m) = %0.2f',x_cell));
        elseif floor(i/Nx) == 3 || i/Nx == 4
            states_title(1,i) = cellstr(sprintf('Froude (-), x(m) = %0.2f',x_cell));
        elseif floor(i/Nx) == 4 || i/Nx == 5
            states_title(1,i) = cellstr(sprintf('Courant Number (-), x(m) = %0.2f',x_cell));
        elseif floor(i/Nx) == 5 || i/Nx == 6
            states_title(1,i) = cellstr(sprintf('Water Surface Elevation (m), x(m) = %0.2f',x_cell));
        end
    end
end
% states_title(1,end+1) = cellstr(sprintf('Water Surface Elevation (m), x(m) = %0.2f',dx*(Nx-1)));
time_string = {'Time (sec)'};
% Table Headers
table_headers = [time_string, states_title];
data_save = zeros(length(time_store),ns*Nx);


if flag_output == 1
    % Detailed Output for each section with all states together
    for i = 1:length(time_store)
        % For all time
        for j = 1:ns
            if j == 2
                ttt = 1;
            end
            % For all states
            for k = 1:Nx
                % For all nodes
                %                 data_table = round(data(i,k,j),decimal_places);
                %                 data_save(i,ns*(k-1) +  j) = data_table;
                data_table = round(data(i,k,j),decimal_places);
                data_save(i,ns*(k-1) +  j) = data_table;
            end
        end
    end
else
    % Detailed Output for each state for each section
    for i = 1:length(time_store)
        % For all time
        for j = 1:ns
            % For all states
            for k = 1:Nx
                % For all nodes
                data_table = round(data(i,k,j),decimal_places);
                data_save(i,k + (j-1)*Nx) = data_table;
            end
        end
    end
end


% Specify the file path including the subfolder
label_plot = strcat(labels.simulation_info.ID,'_',labels.simulation_info.NAME,'_','Detailed_Output.csv');

data_save = [time_save, data_save]; % Concatenating dataset to the time
T = array2table(data_save,'VariableNames',table_headers);
writetable(T,fullfile(folderName,label_plot),'Delimiter',',');
disp('Attention: Data exported in .CSV in Modeling_Results folder');


%% Detailed Output per Cross-Section (Similarly as HEC-RAS)
i_prev = 1;
if flag_elapsed_time == 1
    time_str = 'Elapsed Time (sec)';
else
    time_str = 'Time';
end
Titles_Section = {'x(m)',time_str,' Depth (m)','Discharge (m3/s)','Velocity (m/s)','Froude (-)','Courant Number (-)','WSE (m)'};
data_save_XS = zeros((Nx+1)*length(time_store),6);
for i = 1:Nx
    perc = i/Nx
    % Through each section
    for j = 1:length(time_store)        
        % Through each time
        for k = 1:ns
            row = length(time_store)*(i-1) + j;
            data_save_XS(row,k) = data(j,i,k);

        end
    end
end

zzz = data_save_XS;
clear data_table data_save_XS
for i = 1:Nx
    x_cell = (i-1)*dx;
    if i == 1
        section(1,1) = x_cell;
    end    
    row = length(time_store)*(i-1) + 1;
    row_i = length(time_store)*(i);
    data_save_XS((row + i-1):(row_i + i-1),:) = zzz(row:row_i,:);
    data_save_XS(row_i+1 + i - 1,:) = nan;

    section((row + i-1):(row_i + i-1),:) = x_cell;
    section(row_i+1 + i - 1,:) = NaN;
end
section(size(data_save_XS,1),1) = x_cell;
% section(end:size(data_save_XS,1)) = [];
% section(section == 0) = NaN;
if flag_elapsed_time == 1
    time_vector = time_save;
else
    time_vector = time_begin + time_save/86400; % Days minutes and seconds
end

delta = 0;
for i = 1:Nx
    row = length(time_store)*(i-1) + 1;
    row_i = length(time_store)*(i);
    time_vector_total((row + i-1):(row_i + i-1),1) = time_vector;
    time_vector_total(row_i+1 + i - 1,:) = nan;
end


% Specify the file path including the subfolder
filePath = 'Modeling_Results/Detailed_Output_XS.csv';
label_plot = strcat(labels.simulation_info.ID,'_',labels.simulation_info.NAME,'_','Detailed_Output_XS.csv');

disp('Attention: Data exported in .CSV in Modeling_Results folder');

data_save = [section, time_vector_total, data_save_XS]; % Concatenating dataset to the time
T = array2table(data_save,'VariableNames',Titles_Section);
T.Properties.VariableNames(1:size(data_save,2)) = Titles_Section;
writetable(T,fullfile(folderName,label_plot),'Delimiter',',');
disp('Attention: XS Data exported in .CSV in Modeling_Results folder');
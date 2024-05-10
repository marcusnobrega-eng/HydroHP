%%% Determining Irregular Cross-section Functions %%%
% Developer: Marcus Nobrega Gomes Junior
% Date: 2022/05/03
% Goal - Calculate Hydraulic Properties of Irregular and Regular Sections
% for a given cross-sections and Manning's roughness coefficients

function [y_table, A, P, Rh, y_bar, n_med, Beta, v, B, Q, x_absolute, y,s0] = HP_estimator(flag_plot_HP,dh)
input_table = xlsread('HyProSWE_Input_Data.xlsx','Irregular_Cross_Section');
input_data = input_table(1:5,1);
input_data_coordinates = input_table(2:end,3:end);
flag_length = input_data(1,1); % If == 1, use lengths as main input data, otherwise use absolute values of x (m)
flag_method = input_data(2,1); % If == 1, SCM, else DCM
s0 = input_data(3,1); % Slope in m/m
nm = input_data(4,1); % Main channel roughness
nf = input_data(5,1); % Overbanks channel roughness

if flag_method == 1
    n_channel = input_data_coordinates(1:(end-1),4);
end

% Retrieving Data
x_absolute = input_data_coordinates(:,1);
elevations = input_data_coordinates(:,2);
lengths = input_data_coordinates(1:(end-1),3);
break_point_divider = input_data_coordinates(1:(end),5);

delta = zeros(length(elevations),1);
for i = 1:(length(elevations)-1)
    delta(i) = abs(elevations(i+1,1) - elevations(i,1));
end
delta_h = min(delta(delta > 0));
tic

% Checking input data consistency
if length(elevations) <= 3
    error('Please, enter at least 4 points for elevation and 3 points for manning and lengths. If you have a triangular shape, please enter the invert elevation twice and add a 0 length and 0 manning, such that you have 4 points for elevation and 3 points for manning and lengths')
end

points = (1:1:length(elevations))'; % stations from 1 to n

% Let's assume a maximum 1 cm difference in the depths
% Noise
noise_max = 0.01; % m
% Let's also assume a minimum 0.1 cm difference in the depths, that is, the
% noise
noise_min = 0.001; % m
noise = delta_h/dh; % Noise in m from user input data
if noise > noise_max
    noise = noise_max; % m
elseif noise < noise_min
    noise = noise_min; % m
end

factor = 1; %precision = 1/factor * noise

[au,ia] = unique(elevations,'stable');
Same = ones(size(elevations));
Same(ia) = 0; % repetitive values
noise_i = rand(1,1)*noise;
small_number = noise/100;
% New Elevation and X_values
ii = 0;
for i = 1:(length(elevations) - 1)
    el1 = elevations(i); el2 = elevations(i+1);
    x1 = x_absolute(i); x2 = x_absolute(i+1);
    if el1 == el2 || abs(el1 - el2) == noise
        elevations(i+1) = elevations(i+1) + noise;
        if elevations(i+1) == elevations(i)
            elevations(i+1) = elevations(i+1) + noise;
        end
    end
    if x1 == x2 || abs(x2 - x1) == noise
        x_absolute(i+1) = x_absolute(i+1) + noise;
        if x_absolute(i+1) == x_absolute(i)
            x_absolute(i+1) = x_absolute(i+1) + noise;
        end
    end
end

% if max(isnan(n_channel)) > 0
%     error('Please, enter (n-1) data for Manning coefficient, where n is the number of break-points')
% end

% Roughness Boundary Condition
if flag_method == 1
    n_channel(end+1,1) = 0; % adding last boundary condition
end

% Minimum elevation
min_el = min(elevations); % m
% y (bottom to up)
y = elevations - min_el;
pos_inv = find(y == 0); % position of invert elevation
% If we have more than 1 invert
pos_inv = pos_inv(1);

% x (left to right)
if flag_length == 1
    for i = 1:length(y) % coordinates of each measured point
        if i == 1
            x_absolute(i,1) = 0 + noise;
        else
            x_absolute(i,1) = x_absolute(i-1) + lengths(i-1) + noise;
        end
    end
else % Lengths are already assumed from the input data table
    for i = 1:length(y)
        if i ~= length(y)
            lengths(i) = x_absolute(i+1) - x_absolute(i);
        end
    end
end

%  Alfa min
alfa_min_bound = noise/max(lengths(lengths>1e-8));
big_n = 100000*atan(asin(1)); % big number making sure it is a multiple of 1 rad, so that sin(atan(big_n)) = 1
min_length = min(lengths(lengths>0));

% Invert coordinates
x_invert = x_absolute(pos_inv,1);
y_invert = 0;

% Slopes (taking from x (left-right) y (down-up)
% For point 1 and for the last point
alfa_l = (y(1,1) - y(2,1))/lengths(1,1);

% Unsorted Values
x_left_unsorted = x_absolute(1:(pos_inv-1),1);
y_left_unsorted = y(1:(pos_inv-1),1);
x_right_unsorted = x_absolute(pos_inv + 1:end,1);
y_right_unsorted = y(pos_inv + 1:end,1);
if flag_method == 1
    n_left_unsorted = n_channel(1:(pos_inv-1),1);
    n_right_unsorted = n_channel(pos_inv:(end-1),1);
end

% Maximum depth (left and right)
max_left = max(y_left_unsorted); max_right = max(y_right_unsorted);
max_y = min(max_left, max_right);

% Refreshing values of ymax
pos_r = length(y_right_unsorted);
if max_left ~= max_right
    if max_left > max_y % the maximum is located at left
        z = sort(y_left_unsorted,1,'descend');
        if length(z) == 1 % Case where we have a vertical wall
            z(2,1) = y_invert;
        end
        x_left_first = round(x_absolute(2) - (max_y - z(2))/alfa_l,2);
        % New values of x and y
        x_absolute(1) = x_left_first;
        y(1) = max_y;
        pos_r = length(y_right_unsorted);
    else
        pos_r = find(y_right_unsorted > max_y ,1,'first');
        alfa_r = (y_right_unsorted(pos_r) - y_right_unsorted(pos_r - 1))/lengths(length(y_left_unsorted) + 1 + pos_r-1);
        z = sort(y_right_unsorted,1,'descend');
        x_rigth_last = round(x_absolute(end-1) + (max_y - z(2))/alfa_r,2);
        % New values of x and y
        x_absolute(end) = x_rigth_last;
        y(length(y_left_unsorted) + 1 + pos_r) = max_y;
    end
end
dim = 1:(length(y_left_unsorted) + 1 + pos_r);
y = y(dim,1);
x_absolute = x_absolute(dim,1);
% n_channel = n_channel(dim,1);
points = points(dim);

% New Unsorted Values with New max
x_left_unsorted = x_absolute(1:(pos_inv-1),1);
y_left_unsorted = y(1:(pos_inv-1),1);
x_right_unsorted = x_absolute(pos_inv + 1:end,1);
y_right_unsorted = y(pos_inv + 1:end,1);
if flag_method == 1
    n_left_unsorted = n_channel(1:(pos_inv-1),1);
    n_right_unsorted = n_channel(pos_inv:(end-1),1);
end

% Main Matrix
% table = [points,x_absolute,y,n_channel];

% % Vlookup Function
% Vlookup_eq = @(data,col1,val1,col2) data((find(data(:,col1)==val1,1)),col2); %Vlookup function as Excel
% Vlookup_leq = @(data,col1,val1,col2) data((find(data(:,col1)<=val1,1)),col2); %Vlookup function as Excel

% Sections left
numb_left = length(find(y_left_unsorted >= y_left_unsorted(end)));
% Sections right
numb_right = length(find(y_right_unsorted >= y_right_unsorted(1)));
% Tot sections
tot_sections = numb_left + numb_right - 1; % take one out because both sides are equal

y_l_prev = y_left_unsorted(2:length(y_left_unsorted));
y_l_next = y_left_unsorted(1:(length(y_left_unsorted)-1));

%%%% Precision
precision = 1/factor*noise; % m

%%%% small number >= 1 < 1e-8 + 1
sm = (1e-8 + 1);

%%%% Total_Noise
tot_noise = noise*sum(Same);
% Main loop
i = 0; int_n_p = 0; % integral of n*perimeter

%% Define Main Channel and Overbanks
pos_break = find(break_point_divider == 1); % Position where the divider occurs
% Main Channel Height
ym = y(pos_break); % Main channel height (m)
if pos_break > pos_inv % Left intersection
    % Left intersection
    posm_left = find(y_left_unsorted >= ym,1,'last');
    ym_left_up = y_left_unsorted(posm_left);
    xm_left_up = x_left_unsorted(posm_left);
    ym_left_down = y_left_unsorted(min(posm_left+1,length(y_left_unsorted)));
    xm_left_down = x_left_unsorted(min(posm_left+1,length(y_left_unsorted)));
    % Angles
    if (ym_left_up - ym_left_down <= length(y_left_unsorted)*noise)
        alfa_m_l = big_n;
    else
        alfa_m_l = (ym_left_up - ym_left_down )/(xm_left_down - xm_left_up); % Slope
    end
    xm_left = xm_left_down - (ym - ym_left_down )/alfa_m_l;
    ym_left = ym;
    % Polygons (left - inv - right)
    x_pol = [xm_left; x_left_unsorted((posm_left + 1:end),1); x_invert; x_right_unsorted(1:(pos_break-pos_inv),1)];
    y_pol = [ym_left; y_left_unsorted((posm_left + 1:end),1); y_invert; y_right_unsorted(1:(pos_break-pos_inv),1)];
    % Top-Width
    bm = abs(x_pol(1) - x_pol(end));
    % Area
    am = polyarea(x_pol,y_pol);
    % Perimeter
    polyin = polyshape(x_pol,y_pol);
    pm = perimeter(polyin) - bm; % Taking away the top width
else
    % Right Intersection
    posm_right = find(y_right_unsorted >= ym,1,'first');
    ym_right_up = y_right_unsorted(posm_right);
    xm_right_up = x_right_unsorted(posm_right);
    ym_right_down = y_right_unsorted(max(posm_right-1,1));
    xm_right_down = x_right_unsorted(max(posm_right-1,1));
    % Angles
    if (ym_right_up - ym_right_down < noise*length(y_right_unsorted)) % No depth
        alfa_m_r = big_n;
    else
        alfa_m_r = (ym_right_up - ym_right_down )/(xm_right_up - xm_right_down); % Slope
    end
    xm_right = xm_right_down + (ym - ym_right_down )/alfa_m_r;
    ym_right = ym;
    % Polygons (left - inv - right)
    x_pol = [x_left_unsorted(pos_break:end,1); x_invert; x_right_unsorted(1:(posm_right - 1),1); xm_right];
    y_pol = [y_left_unsorted(pos_break:end,1); y_invert; y_right_unsorted(1:(posm_right - 1),1); ym_right];
    % Top-Width
    bm = abs(x_pol(1) - x_pol(end));
    % Area
    am = polyarea(x_pol,y_pol);
    % Perimeter
    polyin = polyshape(x_pol,y_pol);
    pm = perimeter(polyin) - bm; % Taking away the top width
end
if flag_method ~= 1
    % Number of floodplains
    if pos_break == 1 || pos_break == length(y)
        n_fp = 1;
    else
        n_fp = 2;
    end
end
while i < big_n
    %% Case where i == 1
    i = i + 1;
    n_P_left = 0;
    n_P_right = 0;
    n_P_left_extra = 0;
    n_P_right_extra = 0;
    B_extra = 0;
    P_extra = 0;
    P_extra_left = 0;
    P_extra_right = 0;
    if i == 1 % We are talking about the first point

        %%% Initializing variables
        y_table = 0; h = 0; B = 0; A = 0; Rh = 0; P = 0; Phi = 0; K_c = 0;
        % Look to both sides from pos_inv (invert point)

        % Left Direction
        pos_left = find(y_left_unsorted>sm*y_invert,1,'last');
        y_left_point = y_left_unsorted(pos_left,1);
        x_left_point = x_left_unsorted(pos_left,1);
        if flag_method == 1
            n_left_segment = n_left_unsorted(pos_left,1);
        else
            n_left_segment = nm; % Main channel
        end

        % Right Direction
        pos_right = find(y_right_unsorted>sm*y_invert,1,'first');
        y_right_point = y_right_unsorted(pos_right,1);
        x_right_point = x_right_unsorted(pos_right,1);
        if flag_method == 1
            n_right_segment = n_right_unsorted(pos_right,1);
        else
            n_right_segment = nm; % Main channel
        end

        %%%%%%%%%%%% Angles Calculations %%%%%%%%%%%%
        %%%% Alfa Left %%%%
        % Case 01 - Vertical Point
        if (x_invert - x_left_point <= tot_noise) && (y_left_point - y_invert > tot_noise)
            alfa_l = big_n;
            alfa_l_tang = big_n;
        end
        % Case 02 - Horizontal Point
        if (x_invert - x_left_point > tot_noise) && (y_left_point - y_invert <= tot_noise)
            alfa_l = big_n;
            alfa_l_tang = big_n;
        end
        % Case 03 - Horizontal and Vertical Point
        if (x_invert - x_left_point <= tot_noise) && (y_left_point - y_invert <= tot_noise)
            alfa_l = big_n;
            alfa_l_tang = big_n;
        end
        % Case 04 - Poit with normal slopes
        if (x_invert - x_left_point > tot_noise) && (y_left_point - y_invert > tot_noise)
            alfa_l = (y_left_point - y_invert)/(x_invert - x_left_point);
            alfa_l_tang = alfa_l;
        end

        %%%% Alfa Right %%%%
        % Case 01 - Vertical Point
        if (x_right_point - x_invert <= tot_noise) && (y_right_point- y_invert > tot_noise)
            alfa_r = big_n;
            alfa_r_tang = big_n;
        end
        % Case 02 - Horizontal Point
        if (x_right_point - x_invert > tot_noise) && (y_right_point - y_invert <= tot_noise)
            alfa_r = big_n;
            alfa_r_tang = big_n;
        end
        % Case 03 - Horizontal and Vertical Point
        if (x_right_point - x_invert <= tot_noise) && (y_right_point - y_invert <= tot_noise)
            alfa_r = big_n;
            alfa_r_tang = big_n;
        end
        % Case 04 - Poit with normal slopes
        if (x_right_point - x_invert > tot_noise) && (y_right_point - y_invert > tot_noise)
            alfa_r = (y_right_point - y_invert)/(x_right_point - x_invert);
            alfa_r_tang = alfa_r;
        end

        % Min Angle
        if alfa_l <= alfa_min_bound
            alfa_l_tang = big_n;
        end
        if alfa_r <= alfa_min_bound
            alfa_r_tang = big_n;
        end

        if y_left_point <= y_right_point
            y_moving = y_left_point;
            xleft_point = x_absolute(pos_inv - 1,1);
            precision_section = min(y_left_point - y_invert,precision);
            n_points = floor((y_left_point - y_invert)/(precision_section)); % number of interpolated points
            if n_points == 1 % only one point means no slope
                if x_invert - x_left_point >= sm*noise && alfa_l == big_n
                    P_extra_left = sqrt((x_invert - x_left_point)^2 + (y_invert - y_left_point)^2);
                    n_P_left_extra = P_extra_left*n_left_segment^(3/2);
                    B_extra = (x_invert - x_left_point);
                else
                    B_extra = 0;
                    n_P_left_extra = 0;
                    P_extra_left;
                end
            end
            if n_points == 1 % only one point means no slope
                if x_right_point - x_invert > 1.0001*noise && alfa_r == big_n
                    P_extra_right = sqrt((x_invert - x_right_point)^2 + (y_invert - y_right_point)^2) + B_extra;
                    B_extra = B_extra + (x_right_point - x_invert);
                    n_P_right_extra = (P_extra_right)*n_right_segment^(3/2);
                else
                    n_P_left_extra = 0;
                    P_extra_right = 0;
                end
            end
            P_extra = P_extra_right + P_extra_left;

            %%%%%%%%%%%% Main loop for i == 1 %%%%%%%%%%%%
            for j = 1:(n_points)
                h = precision_section;
                y_table(j+1,1) = y_table(j,1) + h;
                B(j+1,1) = h/alfa_l_tang + h/alfa_r_tang + B(j,1);
                A(j+1,1) = (B(j+1,1) + B(j,1))*h/2 + A(j,1); % Trapezoid
                P(j+1,1) = h/sin(atan(alfa_l_tang)) + h/sin(atan(alfa_r_tang)) + P(j,1);
                Rh(j+1,1) = A(j+1,1)/P(j+1,1);
                Phi(j+1,1) = A(j+1,1)*Rh(j+1,1)^(2/3);
                int_n_p = n_P_left_extra + n_P_right_extra + n_left_segment^(3/2)*h/sin(atan(alfa_l_tang)) + n_right_segment^(3/2)*h/sin(atan(alfa_r_tang)) + int_n_p;
                % Representative Roughness Coefficient
                if flag_method == 1
                    n_med(j+1,1) = (int_n_p/P(j+1,1))^(2/3);
                else
                    if y_table(j+1,1) > ym
                        yf = max(y_table(j+1,1) - ym,0); % Overbank depth
                        af = max(A(j+1,1) - (am + bm*yf),0); % Overbank flow area
                        pf = max(P(j+1,1) - pm,0); % Floodplain perimeter (m)
                        pm_star = max(pm + n_fp*yf,0);
                        am_star = max(am + bm*yf,0);
                        n_med(j+1,1) = (Phi(j+1,1))/(1/nf*af*(af/pf)^(2/3) + 1/nm*am_star*(am_star/pm_star)^(2/3));
                    else
                        yf = 0; % Overbank depth
                        af = 0; % Overbank flow area
                        pf = 0; % Floodplain perimeter (m)
                        pm_star = 0;
                        am_star = 0;
                        n_med(j+1,1) = nm;
                    end
                end
                K_c(j+1,1) = 1/n_med(j+1,1)*Phi(j+1,1);

                if j == (n_points) % final point
                    % Final point - make sure you have the exact surveyed point at the end
                    h_ = y_right_point - y_table(j,1);
                    y_table(j+1,1) = y_table(j,1) + h_;
                    B(j+1,1) = h_/alfa_l_tang + h_/alfa_r_tang + B(j,1) + B_extra;
                    A(j+1,1) = (B(j+1,1) + B(j,1))*h/2 + A(j,1); % Trapezoid
                    P(j+1,1) = h_/sin(atan(alfa_l_tang)) + h_/sin(atan(alfa_r_tang)) + P(j,1) + P_extra;
                    Rh(j+1,1) = A(j+1,1)/P(j+1,1);
                    Phi(j+1,1) = A(j+1,1)*Rh(j+1,1)^(2/3);
                    if n_points == 1
                        int_n_p =  n_left_segment^(3/2)*h/sin(atan(alfa_l_tang)) + n_right_segment^(3/2)*h/sin(atan(alfa_r_tang)) + n_P_right_extra + n_P_left_extra;
                    else
                        int_n_p =  n_left_segment^(3/2)*h/sin(atan(alfa_l_tang)) + n_right_segment^(3/2)*h/sin(atan(alfa_r_tang)) + int_n_p;
                    end
                    % Representative Roughness Coefficient
                    if flag_method == 1
                        n_med(j+1,1) = round((int_n_p/P(j+1,1))^(2/3),3);
                    else
                        if y_table(j+1,1) > ym
                            yf = max(y_table(j+1,1) - ym,0); % Overbank depth
                            af = max(A(j+1,1) - (am + bm*yf),0); % Overbank flow area
                            pf = max(P(j+1,1) - pm,0); % Floodplain perimeter (m)
                            pm_star = max(pm + n_fp*yf,0);
                            am_star = max(am + bm*yf,0);
                            n_med(j+1,1) = round((Phi(j+1,1))/(1/nf*af*(af/pf)^(2/3) + 1/nm*am_star*(am_star/pm_star)^(2/3)),3);
                        else
                            yf = 0; % Overbank depth
                            af = 0; % Overbank flow area
                            pf = 0; % Floodplain perimeter (m)
                            pm_star = 0;
                            am_star = 0;
                            n_med(j+1,1) = nm;
                        end
                    end
                    K_c(j+1,1) = 1/n_med(j+1,1)*Phi(j+1,1);
                end
            end
        else
            x_right_point = x_absolute(pos_inv + 1,1);
            precision_section = min(y_right_point - y_invert,precision);
            n_points = floor((y_right_point - y_invert)/(precision_section)); % number of interpolated points
            if n_points == 1 % only one point means no slope
                if x_right_point - x_invert >= sm*noise && alfa_r == big_n % Additional B_extra
                    P_extra = sqrt((x_right_point - x_invert)^2 + (y_right_point - y_invert)^2);
                    B_extra = x_right_point - x_invert;
                    n_P_right_extra = P_extra*n_right_segment^(3/2);
                else
                    B_extra = 0;
                    n_P_right_extra = 0;
                    P_extra = 0;
                end
            end
            y_moving = y_right_point;
            % For loop to calculate functions
            for j = 1:(n_points)
                h = precision_section;
                B(j+1,1) = h/alfa_l_tang + h/alfa_r_tang + B(j,1);
                y_table(j+1,1) = y_table(j,1) + h;
                A(j+1,1) = (B(j+1,1) + B(j,1))*h/2 + A(j,1); % Trapezoid
                P(j+1,1) = h/sin(atan(alfa_l_tang)) + h/sin(atan(alfa_r_tang)) + P(j,1);
                Rh(j+1,1) = A(j+1,1)/P(j+1,1);
                Phi(j+1,1) = A(j+1,1)*Rh(j+1,1)^(2/3);
                int_n_p = n_P_left_extra + n_P_right_extra + n_left_segment^(3/2)*h/sin(atan(alfa_l_tang)) + n_right_segment^(3/2)*h/sin(atan(alfa_r_tang)) + int_n_p;
                % Representative Roughness Coefficient
                if flag_method == 1
                    n_med(j+1,1) = round((int_n_p/P(j+1,1))^(2/3),3);
                else
                    if y_table(j+1,1) > ym
                        yf = max(y_table(j+1,1) - ym,0); % Overbank depth
                        af = max(A(j+1,1) - (am + bm*yf),0); % Overbank flow area
                        pf = max(P(j+1,1) - pm,0); % Floodplain perimeter (m)
                        pm_star = max(pm + n_fp*yf,0);
                        am_star = max(am + bm*yf,0);
                        n_med(j+1,1) = (Phi(j+1,1))/(1/nf*af*(af/pf)^(2/3) + 1/nm*am_star*(am_star/pm_star)^(2/3));
                    else
                        yf = 0; % Overbank depth
                        af = 0; % Overbank flow area
                        pf = 0; % Floodplain perimeter (m)
                        pm_star = 0;
                        am_star = 0;
                        n_med(j+1,1) = nm;
                    end
                end
                K_c(j+1,1) = 1/n_med(j+1,1)*Phi(j+1,1);
                if j == (n_points) % final point
                    % Final point - make sure you have the exact surveyed point at the end
                    h_ = y_right_point - y_table(j,1);
                    y_table(j+1,1) = y_table(j,1) + h_;
                    B(j+1,1) = h_/alfa_l_tang + h_/alfa_r_tang + B(j,1) + B_extra;
                    A(j+1,1) = (B(j+1,1) + B(j,1))*h/2 + A(j,1); % Trapezoid
                    P(j+1,1) = h_/sin(atan(alfa_l_tang)) + h_/sin(atan(alfa_r_tang)) + P(j,1) + P_extra;
                    Rh(j+1,1) = A(j+1,1)/P(j+1,1);
                    Phi(j+1,1) = A(j+1,1)*Rh(j+1,1)^(2/3);
                    if n_points == 1
                        int_n_p =  n_left_segment^(3/2)*h/sin(atan(alfa_l_tang)) + n_right_segment^(3/2)*h/sin(atan(alfa_r_tang)) + n_P_right_extra + n_P_left_extra;
                    else
                        int_n_p =  n_left_segment^(3/2)*h/sin(atan(alfa_l_tang)) + n_right_segment^(3/2)*h/sin(atan(alfa_r_tang)) + int_n_p;
                    end
                    % Representative Roughness Coefficient
                    if flag_method == 1
                        n_med(j+1,1) = (int_n_p/P(j+1,1))^(2/3);
                    else
                        if y_table(j+1,1) > ym
                            yf = max(y_table(j+1,1) - ym,0); % Overbank depth
                            af = max(A(j+1,1) - (am + bm*yf),0); % Overbank flow area
                            pf = max(P(j+1,1) - pm,0); % Floodplain perimeter (m)
                            pm_star = max(pm + n_fp*yf,0);
                            am_star = max(am + bm*yf,0);
                            n_med(j+1,1) = (Phi(j+1,1))/(1/nf*af*(af/pf)^(2/3) + 1/nm*am_star*(am_star/pm_star)^(2/3));
                        else
                            yf = 0; % Overbank depth
                            af = 0; % Overbank flow area
                            pf = 0; % Floodplain perimeter (m)
                            pm_star = 0;
                            am_star = 0;
                            n_med(j+1,1) = nm;
                        end
                    end
                    K_c(j+1,1) = 1/n_med(j+1,1)*Phi(j+1,1);
                end
            end
        end
        % Previous Positions
        pos_left_previous = pos_left;
        pos_right_previous = pos_right;
    else
        %% Case where i ~= 1

        % Look to left sides from x_point_left and from right side of
        % x_point_right
        y_moving = y_table(end,1); % actual water depth

        % Left Direction
        pos_left = find(y_left_unsorted>sm*y_moving,1,'last');
        y_left_point = y_left_unsorted(pos_left,1);
        x_left_point = x_left_unsorted(pos_left,1);

        % Right Direction
        pos_right = find(y_right_unsorted>sm*y_moving,1,'first');
        y_right_point = y_right_unsorted(pos_right,1);
        x_right_point = x_right_unsorted(pos_right,1);

        % Roughness
        if y_moving <= ym % Inside of the channel
            if flag_method == 1
                n_left_segment = n_left_unsorted(pos_left,1);
                n_right_segment = n_right_unsorted(pos_right,1);
            else
                if (abs(y_left_unsorted(pos_left) - ym) <= noise*length(y_left_unsorted))
                    n_left_segment = nf; % Attention here
                else
                    n_left_segment = nm; % Attention here
                end
                if (abs(y_right_unsorted(pos_right) - ym) <= noise*length(y_right_unsorted))
                    n_right_segment = nf;  % Attention here
                else
                    n_right_segment = nm;  % Attention here
                end
            end
        else % Overbanks
            if flag_method == 1
                n_left_segment = n_left_unsorted(pos_left,1);
                n_right_segment = n_right_unsorted(pos_right,1);
            elseif y_left_unsorted(pos_left) - ym < noise*length(y_left_unsorted)% Check Noises
                n_left_segment = nm; % Attention here
                n_right_segment = nm;  % Attention here
            else
                n_left_segment = nf; % Attention here
                n_right_segment = nf;  % Attention here
            end
        end


        % Checking Discontinuities
        %%% Initializing Varaibles
        Delta_Area_left = 0; Delta_Area_right = 0;
        Delta_B_left = 0; Delta_B_right = 0;
        Delta_P_left = 0; Delta_P_right = 0;

        %%%%%%%%%%%% Angles Calculation %%%%%%%%%%%%
        if pos_left + 1 > length(y_left_unsorted)
            x_prev_left = x_invert;
            y_prev_left = y_invert;
        else
            x_prev_left = (x_left_unsorted(pos_left + 1,1));
            y_prev_left = (y_left_unsorted(pos_left + 1,1));
        end

        %%%% Alfa Left %%%%
        % Case 01 - Vertical Point
        if (x_prev_left - x_left_point <= tot_noise) && (y_left_point - y_prev_left > tot_noise)
            alfa_l = big_n;
            alfa_l_tang = big_n;
        end
        % Case 02 - Horizontal Point
        if (x_prev_left - x_left_point > tot_noise) && (y_left_point - y_prev_left <= tot_noise)
            alfa_l = big_n;
            alfa_l_tang = big_n;
        end
        % Case 03 - Horizontal and Vertical Point
        if (x_prev_left - x_left_point <= tot_noise) && (y_left_point - y_prev_left <= tot_noise)
            alfa_l = big_n;
            alfa_l_tang = big_n;
        end
        % Case 04 - Poit with normal slopes
        if (x_prev_left - x_left_point > tot_noise) && (y_left_point - y_prev_left > tot_noise)
            alfa_l = (y_left_point - y_prev_left)/(x_prev_left - x_left_point);
            alfa_l_tang = alfa_l;
        end
        if pos_right == 1
            x_prev_right = x_invert;
            y_prev_right = y_invert;
        else
            x_prev_right = x_right_unsorted(pos_right - 1,1);
            y_prev_right = y_right_unsorted(pos_right - 1,1);
        end
        %%%% Alfa Right %%%%
        % Case 01 - Vertical Point
        if (x_right_point - x_prev_right <= tot_noise) && (y_right_point- y_prev_right > tot_noise)
            alfa_r = big_n;
            alfa_r_tang = big_n;
        end
        % Case 02 - Horizontal Point
        if (x_right_point - x_prev_right > tot_noise) && (y_right_point - y_prev_right <= tot_noise)
            alfa_r = big_n;
            alfa_r_tang = big_n;
        end
        % Case 03 - Horizontal and Vertical Point
        if (x_right_point - x_prev_right <= tot_noise) && (y_right_point - y_prev_right <= tot_noise)
            alfa_r = big_n;
            alfa_r_tang = big_n;
        end
        % Case 04 - Poit with normal slopes
        if (x_right_point - x_prev_right > tot_noise) && (y_right_point - y_prev_right > tot_noise)
            alfa_r = (y_right_point - y_prev_right)/(x_right_point - x_prev_right);
            alfa_r_tang = alfa_r;
        end

        % Min Angle
        if alfa_l <= alfa_min_bound
            alfa_l_tang = big_n;
        end
        if alfa_r <= alfa_min_bound
            alfa_r_tang = big_n;
        end


        if (pos_left_previous - pos_left) > 1 % More than one movement

            % intersect
            if alfa_l_tang == 0
                x_intersect = x_left_unsorted(pos_left + 1,1);
            else
                x_intersect = x_left_unsorted(pos_left + 1,1) - (y_moving - y_left_unsorted(pos_left + 1,1))/alfa_l;
            end
            x_pol = []; y_pol = [];
            for nn = 1:(pos_left_previous - pos_left)
                x_pol = [x_pol; x_left_unsorted(pos_left_previous - nn + 1)];
                y_pol = [y_pol; y_left_unsorted(pos_left_previous - nn + 1)];
            end
            % Adding intersection
            x_pol = [x_pol;x_intersect];
            y_pol = [y_pol;y_moving];
            % Delta B
            Delta_B_left = abs(x_pol(1) - x_pol(end));
            % Delta A
            Delta_Area_left = polyarea(x_pol,y_pol);
            % Delta P
            polyin = polyshape(x_pol,y_pol);
            Delta_P_left = perimeter(polyin) - Delta_B_left; % Taking away top width
            n_P_left = Delta_P_left*n_left_segment^(3/2);
            % Delta Rh left
            % Phi left
            % Conductance Left
        end


        % Checking Discontinuities
        if (pos_right - pos_right_previous) > 1 % More than one movement
            % intersect
            if alfa_r_tang == 0
                x_intersect = x_right_unsorted(pos_right - 1,1);
            else
                x_intersect = x_right_unsorted(pos_right - 1,1) + (y_moving - y_right_unsorted(pos_right - 1,1))/alfa_r;
            end
            x_pol = []; y_pol = [];
            for nn = 1:(pos_right - pos_right_previous)
                x_pol = [x_pol; x_right_unsorted(pos_right_previous + nn - 1)];
                y_pol = [y_pol; y_right_unsorted(pos_right_previous + nn - 1)];
            end
            % Adding intersection
            x_pol = [x_pol;x_intersect];
            y_pol = [y_pol;y_moving];
            % Delta B
            Delta_B_right = abs(x_pol(1) - x_pol(end));
            % Delta A
            Delta_Area_right = polyarea(x_pol,y_pol);
            % Delta P
            polyin = polyshape(x_pol,y_pol);
            Delta_P_right = perimeter(polyin) - Delta_B_right; % Taking away top width
            % Manning * Perimeter
            n_P_right = Delta_P_right*n_right_segment^(3/2);
        end
        y_moving_end = min(y_right_point,y_left_point);
        %         if (y_moving_end - y_moving)/(precision/100) < 1
        %             error('Please, increase precision. Instability!')
        %         end
        precision_section = min(y_moving_end - y_moving,precision); % meters
        if y_moving_end - y_moving < precision
            ttt = 1;
        end
        n_points = floor((y_moving_end - y_moving)/(precision_section)); % number of interpolated points
        % For loop to calculate functions
        if n_points == 1 % only one point means no slope
            if y_moving_end == y_right_point && y_moving_end == y_left_point && alfa_l == big_n && alfa_r == big_n
                B_extra = x_right_point - x_prev_right + x_prev_left - x_left_point;
                P_extra_left = sqrt((x_prev_left - x_left_point)^2 + (y_prev_left - y_left_point)^2);
                P_extra_right = sqrt((x_right_point - x_prev_right)^2 + (y_right_point - y_prev_right)^2);
            elseif y_moving_end == y_right_point && alfa_r == big_n
                if pos_right == 1
                    P_extra_right = sqrt((x_right_point - x_invert)^2 + (y_right_point - y_invert)^2);
                    B_extra = x_right_point - x_invert;
                else
                    P_extra_right = sqrt((x_right_point - x_prev_right)^2 + (y_right_point -  y_prev_right)^2);
                    B_extra = x_right_point - x_prev_right;
                end
            else % y_moving == y_left
                if pos_left + 1 > length(x_left_unsorted) && alfa_l == big_n
                    P_extra_left = sqrt((x_invert - x_left_point)^2 + (y_invert - y_left_point)^2);
                    B_extra = x_invert - x_left_point;
                elseif alfa_l == big_n
                    P_extra_left = sqrt((x_prev_left - x_left_point)^2 + (y_prev_left - y_left_point)^2);
                    B_extra = x_prev_left - x_left_point;
                end
                % Right
                if pos_right == 1 && alfa_r == big_n
                    P_extra_right = sqrt((x_invert - x_right_point)^2 + (y_invert - y_right_point)^2);
                    B_extra = x_right_point - x_invert + B_extra;
                elseif alfa_r == big_n
                    P_extra_left = sqrt((x_prev_right - x_right_point)^2 + (y_right_point -  y_prev_right^2));
                    B_extra = x_right_point - x_prev_right + B_extra;
                end

            end
            P_extra = P_extra_left + P_extra_right;
            n_P_right_extra = P_extra_right*n_right_segment^(3/2);
            n_P_left_extra = P_extra_left*n_left_segment^(3/2);
        else
            B_extra = 0;
            n_P_right_extra = 0;
            n_P_left_extra = 0;
            P_extra = 0;
            P_extra_left = 0;
            P_extra_right = 0;
        end

        dim_table = length(y_table);
        %%%%%%%%%%%% Main loop for i ~= 1 %%%%%%%%%%%%

        for j = 1:(n_points)
            k = dim_table + j;
            if  j == 1 % We have to add values from discontinuity (Deltas)
                h = precision_section; % meters
                y_table(k,1) = y_table(k-1,1) + h;
                % Roughness
                if y_table(k,1) <= ym % Inside of the channel
                    if flag_method == 1
                        n_left_segment = n_left_unsorted(pos_left,1);
                        n_right_segment = n_right_unsorted(pos_right,1);
                    else
                        if (abs(y_left_unsorted(pos_left) - ym) <= noise*length(y_left_unsorted))
                            n_left_segment = nf; % Attention here
                        else
                            n_left_segment = nm; % Attention here
                        end
                        if (abs(y_right_unsorted(pos_right) - ym) <= noise*length(y_right_unsorted))
                            n_right_segment = nf;  % Attention here
                        else
                            n_right_segment = nm;  % Attention here
                        end
                    end
                else % Overbanks
                    if flag_method == 1
                        n_left_segment = n_left_unsorted(pos_left,1);
                        n_right_segment = n_right_unsorted(pos_right,1);
                    elseif y_left_unsorted(pos_left) - ym < noise*length(y_left_unsorted)% Check Noises
                        n_left_segment = nm; % Attention here
                        n_right_segment = nm;  % Attention here
                    else
                        n_left_segment = nf; % Attention here
                        n_right_segment = nf;  % Attention here
                    end
                end
                B(k,1) = B(k-1,1) + Delta_B_left + Delta_B_right +  h/alfa_l_tang + h/alfa_r_tang;
                A(k,1) = A(k-1,1) + (B(k,1) + B(k-1,1))*h/2 + Delta_Area_left + Delta_Area_right;
                P(k,1) = h/sin(atan(alfa_l_tang)) + h/sin(atan(alfa_r_tang)) + P(k-1,1) + Delta_P_left + Delta_P_right;
                Rh(k,1) = A(k,1)/P(k,1);
                Phi(k,1) = A(k,1)*Rh(k,1)^(2/3);
                int_n_p =  n_P_left + n_P_right + n_P_right_extra + n_P_left_extra + n_left_segment^(3/2)*h/sin(atan(alfa_l_tang)) + n_right_segment^(3/2)*h/sin(atan(alfa_r_tang)) + int_n_p;
                % Representative Roughness Coefficient
                if flag_method == 1
                    n_med(k,1) = (int_n_p/P(k,1))^(2/3);
                else
                    if y_table(k,1) > ym
                        yf = max(y_table(k,1) - ym,0); % Overbank depth
                        af = max(A(k,1) - (am + bm*yf),0); % Overbank flow area
                        pf = max(P(k,1) - pm,0); % Floodplain perimeter (m)
                        pm_star = max(pm + n_fp*yf,0);
                        am_star = max(am + bm*yf,0);
                        n_med(k,1) = (Phi(k,1))/(1/nf*af*(af/pf)^(2/3) + 1/nm*am_star*(am_star/pm_star)^(2/3));
                    else
                        yf = 0; % Overbank depth
                        af = 0; % Overbank flow area
                        pf = 0; % Floodplain perimeter (m)
                        pm_star = 0;
                        am_star = 0;
                        n_med(k,1) = nm;
                    end
                end
                K_c(k,1) = 1/n_med(k,1)*Phi(k,1);
            else
                % Functions in terms of depth
                h = precision_section;
                y_table(k,1) = h + y_table(k-1,1);
                % Roughness
                if y_table(k,1) <= ym % Inside of the channel
                    if flag_method == 1
                        n_left_segment = n_left_unsorted(pos_left,1);
                        n_right_segment = n_right_unsorted(pos_right,1);
                    else
                        if (abs(y_left_unsorted(pos_left) - ym) <= noise*length(y_left_unsorted))
                            n_left_segment = nf; % Attention here
                        else
                            n_left_segment = nm; % Attention here
                        end
                        if (abs(y_right_unsorted(pos_right) - ym) <= noise*length(y_right_unsorted))
                            n_right_segment = nf;  % Attention here
                        else
                            n_right_segment = nm;  % Attention here
                        end
                    end
                else % Overbanks
                    if flag_method == 1
                        n_left_segment = n_left_unsorted(pos_left,1);
                        n_right_segment = n_right_unsorted(pos_right,1);
                    elseif y_left_unsorted(pos_left) - ym < noise*length(y_left_unsorted)% Check Noises
                        n_left_segment = nm; % Attention here
                        n_right_segment = nm;  % Attention here
                    else
                        n_left_segment = nf; % Attention here
                        n_right_segment = nf;  % Attention here
                    end
                end
                B(k,1) = h/alfa_l_tang + h/alfa_r_tang + B(k-1,1);
                A(k,1) = (B(k,1) + B(k-1,1))*h/2 + A(k-1,1); % Trapezoid
                P(k,1) = h/sin(atan(alfa_l_tang)) + h/sin(atan(alfa_r_tang)) + P(k-1,1);
                Rh(k,1) = A(k,1)/P(k,1);
                Phi(k,1) = A(k,1)*Rh(k,1)^(2/3);
                int_n_p = n_left_segment^(3/2)*h/sin(atan(alfa_l_tang)) + n_right_segment^(3/2)*h/sin(atan(alfa_r_tang)) + int_n_p;
                % Representative Roughness Coefficient
                if flag_method == 1
                    n_med(k,1) = (int_n_p/P(k,1))^(2/3);
                else
                    if y_table(k,1) > ym
                        yf = max(y_table(k,1) - ym,0); % Overbank depth
                        af = max(A(k,1) - (am + bm*yf),0); % Overbank flow area
                        pf = max(P(k,1) - pm,0); % Floodplain perimeter (m)
                        pm_star = max(pm + n_fp*yf,0);
                        am_star = max(am + bm*yf,0);
                        n_med(k,1) = (Phi(k,1))/(1/nf*af*(af/pf)^(2/3) + 1/nm*am_star*(am_star/pm_star)^(2/3));
                    else
                        yf = 0; % Overbank depth
                        af = 0; % Overbank flow area
                        pf = 0; % Floodplain perimeter (m)
                        pm_star = 0;
                        am_star = 0;
                        n_med(k,1) = nm;
                    end
                end
                K_c(k,1) = 1/n_med(k,1)*Phi(k,1);
            end

            if j == (n_points) % final point
                % Final point - make sure you have the exact surveyed point at the end
                h_ = y_moving_end - y_table(k-1,1);
                y_table(k,1) = y_table(k-1,1) + h_;
                % Roughness
                if y_table(k,1) <= ym % Inside of the channel
                    if flag_method == 1
                        n_left_segment = n_left_unsorted(pos_left,1);
                        n_right_segment = n_right_unsorted(pos_right,1);
                    else
                        if (abs(y_left_unsorted(pos_left) - ym) <= noise*length(y_left_unsorted))
                            n_left_segment = nf; % Attention here
                        else
                            n_left_segment = nm; % Attention here
                        end
                        if (abs(y_right_unsorted(pos_right) - ym) <= noise*length(y_right_unsorted))
                            n_right_segment = nf;  % Attention here
                        else
                            n_right_segment = nm;  % Attention here
                        end
                    end
                else % Overbanks
                    if flag_method == 1
                        n_left_segment = n_left_unsorted(pos_left,1);
                        n_right_segment = n_right_unsorted(pos_right,1);
                    elseif y_left_unsorted(pos_left) - ym < noise*length(y_left_unsorted)% Check Noises
                        n_left_segment = nm; % Attention here
                        n_right_segment = nm;  % Attention here
                    else
                        n_left_segment = nf; % Attention here
                        n_right_segment = nf;  % Attention here
                    end
                end
                B(k,1) = h_/alfa_l_tang + h_/alfa_r_tang + B(k-1,1) + B_extra;
                A(k,1) = (B(k,1) + B(k-1,1))*h_/2 + A(k-1,1); % Trapezoid
                P(k,1) = h_/sin(atan(alfa_l_tang)) + h_/sin(atan(alfa_r_tang)) + P(k-1,1) + P_extra;
                Rh(k,1) = A(k,1)/P(k,1);
                Phi(k,1) = A(k,1)*Rh(k,1)^(2/3);
                int_n_p = n_left_segment^(3/2)*h/sin(atan(alfa_l_tang)) + n_right_segment^(3/2)*h/sin(atan(alfa_r_tang)) + int_n_p;
                % Representative Roughness Coefficient
                if flag_method == 1
                    n_med(k,1) = (int_n_p/P(k,1))^(2/3);
                else
                    if y_table(k,1) >= ym
                        yf = max(y_table(k,1) - ym,0); % Overbank depth
                        af = max(A(k,1) - (am + bm*yf),0); % Overbank flow area
                        pf = max(P(k,1) - pm,0); % Floodplain perimeter (m)
                        pm_star = max(pm + n_fp*yf,0);
                        am_star = max(am + bm*yf,0);
                        n_med(k,1) = (Phi(k,1))/(1/nf*af*(af/pf)^(2/3) + 1/nm*am_star*(am_star/pm_star)^(2/3));
                    else
                        yf = 0; % Overbank depth
                        af = 0; % Overbank flow area
                        pf = 0; % Floodplain perimeter (m)
                        pm_star = 0;
                        am_star = 0;
                        n_med(k,1) = nm;
                    end
                end
                K_c(k,1) = 1/n_med(k,1)*Phi(k,1);
            end
        end
        % Previous Positions
        pos_left_previous = pos_left;
        pos_right_previous = pos_right;
    end
    % Checking i
    if round(y_table(end),3) == round(max_y,3) % Stop de algorithm
        i = big_n;
    end
end

% Centroid Coordinates
int_a_y = 0; % Integral of A(y)dy
for i = 1:(length(A))
    if i == 1
        y_bar(i,1) = 0;
        int_a_y(i,1) = 0;
    else
        int_a_y(i,1) = (A(i) - A(i-1))*(y_table(i) + y_table(i-1))/2 + int_a_y(i-1);
        y_bar(i,1) = int_a_y(i,1)/A(i,1);
    end
end

% Flow Discharge Calculations
Q = K_c*sqrt(s0);

% Velocity
v = Q./A; % m/s

% Beta - Boussinesq factor
kappa = 0.41;
g = 9.81; % m/s2
Beta = (1 + (g*n_med.^2)./(Rh.^(1/3)*kappa^2));

%% Plotting Results
% Plotting Channel
if flag_plot_HP == 1
close all
subplot(1,2,1)
set(gcf,'units','inches','position',[4,4,6.5,4])
mark_size = 5;
plot(x_absolute,y,'linewidth',2,'color','black')
xlabel('x ($m$)','Interpreter','latex');
ylabel('y ($m$)','Interpreter','latex');
xlim([min(x_absolute) max(x_absolute)])
grid on
hold on
scatter(x_absolute,y,'black')
subplot(1,2,2)
n_med(1,1) = inf;
plot(n_med(2:end,1),y_table(2:end,1),'linewidth',2,'color','black')
xlabel('Manning`s coefficient (SI)','Interpreter','latex');
ylabel('y ($m$)','Interpreter','latex');
xlim([0.9*min(n_med) 1.1*max(n_med(~isinf(n_med)))])
grid on
exportgraphics(gcf,'Cross_Section.pdf','ContentType','vector')


subplot(2,4,1)
set(gcf,'units','inches','position',[4,2,7.5,5])
sz = 5;
c = linspace(1,sz,length(y_table));
scatter(A,y_table,sz,c,'filled')
grid on
grid on
xlabel('Area ($m^2$)','Interpreter','latex');
ylabel('y ($m$)','Interpreter','latex');
% xlim([0 4])
subplot(2,4,2)
grid on
scatter(P,y_table,sz,c,'filled')
grid on
xlabel('Perimeter ($m$)','Interpreter','latex');
ylabel('y ($m$)','Interpreter','latex');
% xlim([0 4])
subplot(2,4,3)
grid on
scatter(Rh,y_table,sz,c,'filled')
grid on
xlabel('Hydraulic Radius ($m$)','Interpreter','latex');
ylabel('y ($m$)','Interpreter','latex');
% xlim([0 4])
subplot(2,4,4)
grid on
scatter(B,y_table,sz,c,'filled')
grid on
xlabel('Top width ($m$)','Interpreter','latex');
ylabel('y ($m$)','Interpreter','latex');
subplot(2,4,5)
grid on
scatter(K_c,y_table,sz,c,'filled')
grid on
xlabel('Conveyance ($m^3/s$)','Interpreter','latex');
ylabel('y ($m$)','Interpreter','latex');
subplot(2,4,6)
sz = 5;
c = linspace(1,sz,length(y_table));
scatter(Phi,y_table,sz,c,'filled')
grid on
xlabel('$\Phi$ ($m^{5/3}$)','Interpreter','latex');
ylabel('y ($m$)','Interpreter','latex');
subplot(2,4,7)
scatter(y_bar,y_table,sz,c,'filled')
grid on
xlabel('$\bar{y}$ (m)','Interpreter','latex');
ylabel('y ($m$)','Interpreter','latex');
subplot(2,4,8)
scatter(Q,y_table,sz,c,'filled')
grid on
xlabel('Flow discharge ($m^3/s$)','Interpreter','latex');
ylabel('y ($m$)','Interpreter','latex');
exportgraphics(gcf,'Hydraulic_Properties.pdf','ContentType','vector')
toc

% Rating Curve
close all
subplot(3,1,1)
set(gcf,'units','inches','position',[4,4,6.5,4])
mark_size = 5;
plot(x_absolute,y,'linewidth',2,'color','black')
xlabel('x ($m$)','Interpreter','latex');
ylabel('y ($m$)','Interpreter','latex');
xlim([min(x_absolute) max(x_absolute)])
grid on
subplot(3,1,2)
scatter(Q,y_table,sz,c,'filled')
xlabel('Flow discharge ($m^3/s$)','Interpreter','latex');
ylabel('y ($m$)','Interpreter','latex');
grid on
box on
% Velocity
subplot(3,1,3)
scatter(Q./A,y_table,sz,c,'filled')
xlabel('Velocity ($m/s$)','Interpreter','latex');
ylabel('y ($m$)','Interpreter','latex');
grid on
box on
exportgraphics(gcf,'Rating Curve.pdf','ContentType','vector')

% Plotting Normalized Values
set(gcf,'units','inches','position',[4,2,8,4])
subplot(1,5,1)
scatter(Q/max(Q),y_table/max(y_table),sz,c,'filled')
xlabel('$Q/Q_p$','Interpreter','latex');
ylabel('$y/y_{max} $','Interpreter','latex');
title(['$Q_p (m^3/s) = $ ',num2str(round(max(Q),2))],'interpreter','latex')
axis equal
grid on
xlim([0 1]); ylim([0 1]);
subplot(1,5,2)
scatter(A/max(A),y_table/max(y_table),sz,c,'filled')
xlabel('$A/A_{max}$','Interpreter','latex');
ylabel('$y/y_{max}$','Interpreter','latex');
title(['$A_{max} (m^2) = $ ',num2str(round(max(A),2))],'interpreter','latex')
axis equal
grid on
xlim([0 1]); ylim([0 1]);
subplot(1,5,3)
scatter(Phi/max(Phi),y_table/max(y_table),sz,c,'filled')
xlabel('$\Phi/\Phi_{max}$','Interpreter','latex');
ylabel('$y/y_{max} $','Interpreter','latex');
title(['$\Phi_{max} (m^2) = $ ',num2str(round(max(Phi),2))],'interpreter','latex')
axis equal
grid on
xlim([0 1]); ylim([0 1]);
subplot(1,5,4)
scatter(K_c/max(K_c),y_table/max(y_table),sz,c,'filled')
xlabel('$K_c/K_{c,max}$','Interpreter','latex');
ylabel('$y/y_{max} $','Interpreter','latex');
title(['$K_{c,max} (m^3/s) = $ ',num2str(round(max(K_c),2))],'interpreter','latex')
axis equal
grid on
xlim([0 1]); ylim([0 1]);
subplot(1,5,5)
scatter((Q./A)/(max(Q./A)),y_table/max(y_table),sz,c,'filled')
xlabel('$v/v_{c,max}$','Interpreter','latex');
ylabel('$y/y_{max} $','Interpreter','latex');
title(['$v_{max} (m/s) = $ ',num2str(round(max(Q./A),2))],'interpreter','latex')
axis equal
grid on
xlim([0 1]); ylim([0 1]);
exportgraphics(gcf,'Normalized_Values.pdf','ContentType','vector')
close all

end
end

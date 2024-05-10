%%% --------- HydroHP Model ----------- %%%
% Developer: Marcus Nobrega Gomes Junior
% 5/1/2023
% Goal: Solution of 1-D SVE for given cross-section functions of Area, Perimeter, and
% top Width
% If you have any issues, please contact me at
% marcusnobrega.engcivil@gmail.com


% -------------------- All Rights Reserved ---------------------- %
        
clear all
clc
warning('off') % Deactivate Warnings

% ------------------- Main Input Table ------------------- %
HydroHP_Input_File = 'HydroHP_Input_Data.xlsx'; % Enter the name of the HydroHP-1D input file

%% 1.0 -  Pre-Processing 
% Reading the Input Data
Read_Input_Data % Here we read the .xlsx input data file. Please don't change the name of this file.

% Checking if at least one boundary condition is considered
if flag_hydrograph ~= 1 && flag_nash ~= 1 && flag_stage_hydrograph ~= 1 && flag_outlet ~= 0
    error('Please enter at least 1 internal boundary condition.')
end

% Checking if there is conflicting boundary conditions
if flag_hydrograph == 1 && flag_nash == 1
    error('Please choose either an observed inflow hydrograph entered in a tabular format or a nash-type hydrograph.')
end

% Checking if there is conflicting cross section
if flag_section > 4
    error('Please, enter a the index indicating which type of cross-section is being simulated. Read the instruction in the .xlsx file')
end

% Checking if there is conflicting cross section
if flag_stage_hydrograph == 1 && flag_hydrograph == 1
    error('Please, the inlet can only have either a stage hydrograph or a flow hydrograph')
end
%% 2.0 - Initial Boundary Conditions
% ------------- Inflow Hydrograph ------------- %
if flag_hydrograph == 1
    % We already read the hydrograph in Read_Input_Data file
elseif flag_nash == 1
    % 2nd option - Model the hydrograph using a nash function
    %%% Q(t) = Qb(t) + (Qp(t) - Qb(t))*(t/TP*EXP(1 - t/TP))^Beta
    Inflow_Hydrograph_fun = @(t)(Qb + (Qp - Qb).*(t/(Tp*3600).*exp(1 - (t)/(Tp*3600))).^Beta);
    time = [0 tf]'; % begin and end in min
else
    time = [0 tf]'; % begin and end in min
end

if flag_stage_hydrograph == 1
    Stage_Hydrograph = he1;
end

% -------------  Outlet Boundary Condition ------------- %
% flag_outlet = 1; % 1 = normal depth, flag_outlet >< 1, stage hydrograph
% at the outlet following a wave function

if flag_outlet ~= 1
    %%% Wave Properties for Outlet Stage Hydrograph
    %     x_wave = L_wave/1; % point position in wave x direction;
    k_wave = 2*pi/L_wave;
    sigma_wave = 2*pi./(T_wave*3600);
    h_wave_function = @(t)(h_0_wave + H_0_wave/2.*cos(k_wave.*x_wave - sigma_wave*t));
end

% Time Calculations
time = time*60; % time in seconds
[a1,~] = size(time); % Length of time
tt_h = time(a1,1); % End of hydrograph in seconds
tt = min(tf*60,tt_h); % End of simulation in seconds
Nt = tt/dt; % Number of time-steps in the simulations

% Recording Times
time_records_min = animation_time; % Minutes
time_store = [0:time_records_min*60:tt]; % number of steps necessary to reach the recording vector
Nat = time_records_min*60/dt; % Number of time-steps within an animation time
tint = linspace(0,tt,Nt); % Generate Nt points within 0 and tt(sec)

time_save = zeros(length(time_store),1); % Time
Flow_Area = zeros(length(time_store),Nx); % Flow area
Discharge = zeros(length(time_store),Nx); % Flow discharge
Depth = zeros(length(time_store),Nx); % Depth
Velocity = zeros(length(time_store),Nx); % Velocity
Froude = zeros(length(time_store),Nx); % Froude
Courant = zeros(length(time_store),Nx); % Courant number

if flag_hydrograph == 1
    Qe1int = max(interp1(time,Qe1(:,1),tint,'pchip'),0); % Interpolated flow 
    % Assuming no negative flows
    Qe1int = Qe1int';
elseif flag_nash == 1
    Qe1int = Inflow_Hydrograph_fun(tint)';
else
    tiny_flow = 1e-8;
    Qe1int = tiny_flow*ones(1,length(tint)); % No inflow hydrograph
end

if flag_stage_hydrograph == 1
    he1int = max(interp1(time_stage*60,he1(:,1),tint,'pchip'),0); % Interpolated depth 
    he1int = he1int';
end

%% 3.0 - Pre-Allocation of Arrays

% Channel Discretization
dx = L/(Nx-1); % Channel discretization length in meters

% Friction Data
flag_friction = 1; % If 1, Manning, otherwise DW

% Manning
if flag_manning ~= 1
    nm = repmat(nm,Nx,1); % Bottom slope in m/m for all reaches
else
    nm = manning_values; % Space varying Manning's
end

% Pre-allocating arrays
% Matrices
x = (0:dx:L)'; % x discretization in meters
y = zeros(Nx,2);
q1 = zeros(Nx,2);
q2 = zeros(Nx,2);
f1 = zeros(Nx,2);
f2 = zeros(Nx,2);
J2 = zeros(Nx,2);
q1_back = q1(1:(Nx-2),2);
q1_forward = zeros(Nx-2,2);
q2_back = zeros(Nx-2,2);
q2_forward = zeros(Nx-2,2);
f1_back = zeros(Nx-2,2);
f1_forward = zeros(Nx-2,2);
f2_back = zeros(Nx-2,2);
f2_forward = zeros(Nx-2,2);
J2_back = zeros(Nx-2,2);
J2_forward = zeros(Nx-2,2);
ybar = zeros(Nx,2);
Fr = zeros(Nx,2);
Cn = zeros(Nx,2);

%% 4.0 Channel Data (Cross Section)
% Slope
if flag_slope ~= 1 && flag_elevation ~= 1
    I0 = repmat(I0,(Nx-1),1); % Bottom slope in m/m for all reaches. This is only valid for closed-form sections
elseif flag_slope == 1
    I0 = bottom_slope; % From read input data script
end

if flag_elevation == 1 % We are entering the elevations of each node
    for i = 1:(Nx-1)
        if i+1 > length(inv_el)
            error('Please make sure to add enough invert elevation data.')
        end
        I0(i,1) = (inv_el(i+1) - inv_el(i))/dx;
    end
end

% Outlet Slope
if flag_outlet == 1
    I0(end+1) = s_outlet;
end

% Space-varying trapezoid cross-section data
if flag_trapezoid_sections == 1
    b = cross_section.trapezoid.b;
    Z1 = cross_section.trapezoid.slope_left;
    Z2 = cross_section.trapezoid.slope_right;
end

% Intializing channel data
sm = 1e-12; % Small number
b = sm + b; Z1 = sm + Z1; Z2 = sm + Z2; D = sm + D; a = sm + a;
% flag_section - If 1, trapezoid, if 2, circular, if 3, paraboloid, if 4 - Irregular

% Invert Elevations
if flag_elevation ~=1
    inv_el = zeros(Nx,1);
    for i = 1:Nx
        if i == 1
            inv_el(i) = el;
        else
            inv_el(i) = inv_el(i-1) - (I0(i-1)*dx);
        end
    end
end

% -------------  Geometrical Functions for all Cros-Sections ------------ %
syms b_ y_ Z1_ Z2_ Q_ I0_ D_ a_
dim_all = 1e-6*(y_ + Z1_ + Z2_ + a_ + D_ + b_);
if flag_section == 1
    B = b_ + y_.*(Z1_ + Z2_) + + dim_all; % user defined function (top width)
    B_function = matlabFunction(B);
    P = b_ + y_.*(sqrt(1 + Z1_^2) + sqrt(1 + Z2_^2)) + dim_all; % Perimeter Function % user defined function
    P_function = matlabFunction(P);
    A = (2*b_ + y_.*(Z1_ + Z2_))*y_/2 + dim_all; % Area function % user defined function
    A_function = matlabFunction(A); % Function describing the area in terms of y
    centroid = y_ - int(A,y_)./A + dim_all; % 1st order momentum
    ybar_function = matlabFunction(centroid); % Function describing ybar in terms of y
end
if flag_section == 2
    % Circular Section
    theta = 2*acos(1 - 2.*y_./D_) + dim_all;
    B = D_.*sin(theta/2) ; % top width
    B_function = matlabFunction(B);
    P = theta.*D_/2 ; % perimeter
    P_function = matlabFunction(P);
    A = D_.^2/8.*(theta - sin(theta)) ; % area
    A_function = matlabFunction(A); % Function describing the area in terms of y
    Ybar = y_ - (D_.*(- cos(theta/2)/2 + 2.*sin(theta/2).^3./(3*(theta - sin(theta))))); % Very much attention here
    ybar_function = matlabFunction(Ybar);
end

if flag_section == 3
    % Parabolic Section
    % Area Function
    A = 4.*(y_.^3/2)./(3*sqrt(a_)) + dim_all; % m2
    A_function = matlabFunction(A); % Function describing the area in terms of y
    % Top Width
    B = 3/2.*A./y_ + dim_all; % m
    B_function = matlabFunction(B);
    % Hydraulic Perimeter
    P = dim_all + sqrt(y_)./sqrt(a_).*(sqrt(1 + 4*a_.*y_) + 1./(2*a_).*(log(2*sqrt(a_).*sqrt(y_) + sqrt(1 + 4*a_.*y_))));
    P_function = matlabFunction(P);
    Y_bar = y_ - 2/5*y_ + dim_all;
    ybar_function = matlabFunction(Y_bar);
end

if flag_section ~= 4
    %%%%%%% Hydraulic Radius %%%%%%%
    Rh = A/P; % Hydraulic Radius Function
    Rh_function = matlabFunction(Rh); % Function describing the hydraulic radius in terms of y
end

% Vlookup Function
Vlookup_eq = @(data,col1,val1,col2) data((find(data(:,col1)==val1,1,'first')),col2); %Vlookup function as Excel
Vlookup_l = @(data,col1,val1,col2) data((find(data(:,col1)<val1,1,'last')),col2); %Vlookup function as Excel]
Vlookup_g = @(data,col1,val1,col2) data((find(data(:,col1)>val1,1,'first')),col2); %Vlookup function as Excel
fv = 1 + 1e-4; % Factor to avoid fails in vlookup function
min_area = 0.01;
min_depth = 0.01;

% Minimum Value
if flag_section == 4
    min_depth = 0.02; % m
    min_area = Vlookup_l(irr_table,1,min_depth*fv,2); 
end


% Initial Guess
if flag_section == 1
    y0_guess = 1;
elseif flag_section == 2
    y0_guess = D/2;
elseif flag_section == 3
    y0_guess = 1;
end

%% 5.0- Initial Values for Simulation
Q0 = Qe1int(1,1); % Flow at inlet section at time 0 comes from the table
if flag_stage_hydrograph == 1
    h0 = he1int(1,1); % Water depth at x = 0 at time = 0
end
if flag_friction == 1
    if flag_section ~= 4
        if Q0 == 0
            Q0 = sm; % Numerical Constraint
        end
            y0 = uniformeM(nm,Q0,b,Z1,Z2,a,D,I0,P,A,y0_guess); % normal depth using manning equation
            if imag(y0(1)) ~= 0
                y0 = uniformeM(nm,Q0,b,Z1,Z2,a,D,I0,P,A,5*y0_guess);
                if imag(y0(1)) ~= 0
                    y0 = uniformeM(nm,Q0,b,Z1,Z2,a,D,I0,P,A,10*y0_guess);
                end
            end
        % Stage_Hydrograph Boundary Condition
        if flag_stage_hydrograph == 1
            y0(1,1) = h0;
        end
        % More Initial Boundary Conditions for Area, Velocity, Perimeter and Rh
        A0 = A_function(D,Z1,Z2,a,b,y0); % Cross section area in m2
        u0 = (Q0./A0)'; % Initial velocity in m/s
        P0 = P_function(D,Z1,Z2,a,b,y0); % Hydraulic perimeter in m
        Rh0 = A0./P0; % Hydraulic radius at time 0
        % Boundary Conditions
        y(:,1) = y0; % all sub-reaches with y0 at the beginning
        q1(:,1) = A0; % all sub-reaches with same area A0 at the beginning
        q2(:,1) = Q0; % Assuming permanent conditions at the beginning
        f1(:,1) = q2(:,1);
        % f2 depends on ybar
    else % Irregular Cross-Section
        % [y_irr, A_irr, P_irr, Rh_irr, y_bar_irr, n_med_irr, Beta_irr, u_irr, B_irr, Q_irr];
        % [   1,    2,     3,      4,    5,            6,          7,     8,      9,    10]

        if max(irr_table(:,10)) == 0 % No outflow and S = 0
            % Here we are modeling a channel with no slope
            % We search Everything Using the Depth instead of the Flow
            col1 = 1; % Searching with the Col of Y
            % Stage_Hydrograph Boundary Condition
            if flag_stage_hydrograph == 1
                y0(1,1) = max(h0,irr_table(2,1));
                s_v = y0; % Searching Variable
            else
                error('Please, add a minimum slope value or enter a stage-hydrograph boundary condition.')
            end
            y0 = Vlookup_g(irr_table,col1,s_v*fv,1);
            A0 = Vlookup_g(irr_table,col1,s_v*fv,2);
            P0 = Vlookup_g(irr_table,col1,s_v*fv,3);
            Rh0 = Vlookup_g(irr_table,col1,s_v*fv,4);

        else % Now we are modeling a channel with slope
            if (flag_hydrograph == 1 || flag_nash == 1) && flag_stage_hydrograph ~= 1
                col1 = 10; % Col with Q
                Q_min_table = irr_table(3,10); % ATTENTION HERE
                s_v = max(Q0,Q_min_table); % Searching Variable
            elseif flag_stage_hydrograph == 1
                col1 = 1; % Col with y or h
                h_min_table = irr_table(3,1);
                s_v = max(h_min_table,h0); % Searching Variable
            elseif flag_outlet == 0
                col1 = 1; % Col with y or h
                h_min_table = irr_table(3,1);
                s_v = max(h_min_table,0); % Searching Variable
                Q_min_table = irr_table(3,10);
                Q0 = Q_min_table;
            end
            Q0 = max(irr_table(2,end),Q0); % Allowing minimum value of Q0 larger than 0
            y0 = Vlookup_l(irr_table,col1,s_v*fv,1);
            A0 = Vlookup_l(irr_table,col1,s_v*fv,2);
            P0 = Vlookup_l(irr_table,col1,s_v*fv,3);
            Rh0 = Vlookup_l(irr_table,col1,s_v*fv,4);
        end
        % Boundary Conditions
        y(:,1) = y0; % all sub-reaches with y0 at the beginning
        q1(:,1) = A0; % all sub-reaches with same area A0 at the beginning
        q2(:,1) = Q0; % Assuming permanent conditions at the beginning
        f1(:,1) = q2(:,1);
    end
    if flag_outlet ~=1 % Bay or Ocean Boundary Condition
        % Stage Hydrograph Boundary Condition
        time_wave = 0; % time in seconds
        y(Nx,1) = h_wave_function(time_wave);
        if flag_section ~= 4
            q1(Nx,1) = A_function(D,Z1(Nx),Z2(Nx),a,b(Nx),y(Nx,1));
        else
            % We search Everything Using the Depth instead of the Flow
            col1 = 1; % Searching with the Col of Flow
            q1(Nx,1) = Vlookup_g(irr_table,col1,y(Nx,1)*fv,2);
        end
    end

    % Hydraulic Radius
    if flag_section ~= 4
        Rh_outlet = Rh_function(D,Z1(Nx),Z2(Nx),a,b(Nx),y(Nx,1));
    else
        for mm = 1:(length(irr_table(1,:))-1)
            interp_base = q1(Nx,1); % Value that will be used for interpolation (area)
            area_smaller = Vlookup_l(irr_table,2,interp_base,2); % Smaller values
            if isempty(area_smaller)
                area_smaller = 0;
            end
            area_larger = Vlookup_g(irr_table,2,interp_base,2); % Larger values
            col1 = 2; % Interpolating from area values
            if interp_base <= min_area
                var_outlet(mm,1,1) = irr_table(2,mm); % Smaller values
            else
                var_outlet(mm,1,1) = Vlookup_l(irr_table,col1,interp_base,mm); % Smaller values
            end
            var_outlet(mm,1,2) = Vlookup_g(irr_table,col1,interp_base,mm); % Larger values            
            alfa_var_outlet(mm,1) = sqrt((interp_base - area_smaller)/(area_larger - area_smaller));
        end


        % [y_table, A, P, Rh, y_bar, n_med, Beta, v, B, Q]
        % [   1,    2, 3, 4,    5,    6,      7,  8,  9, 10]
        col_var = 4; % Calculating Hydraulic Radius
        % Var* = Var(-) + alfa*(Var(+) - Var(-))
        Rh_outlet = var_outlet(col_var,1,1) + alfa_var_outlet(col_var,1)*(var_outlet(col_var,1,2) - var_outlet(col_var,1,1)); % Interpolated Hydraulic Radius
        % Var* = Var(-) + alfa*(Var(+) - Var(-))
        col_var = 6;
        nm(end,1) = var_outlet(col_var,1,1) + alfa_var_outlet(col_var,1)*(var_outlet(col_var,1,2) - var_outlet(col_var,1,1));
    end

        if flag_outlet == 1
            u = (1./nm(Nx)).*Rh_outlet^(2/3)*I0(Nx)^0.5; % Normal depth at the outlet
            flow_dir = 1;
        else
            wse_dif = y(Nx-1,2-1) + inv_el(Nx-1) - y(Nx,2) - inv_el(Nx); % Difference in wse
            out_slope = abs(wse_dif)/dx; % Friction slope at the outlet as a diffusive model
            if wse_dif < 0 
                ttt = 1;
            end
              
            if flag_stage_hydrograph ~= 1 && flag_nash ~= 1 && flag_hydrograph ~= 1
                % Only Outlet Tidal B.C.
                if wse_dif > 0 && y(Nx-1,2-1) <= fv*1e-3
                    out_slope = 0;
                end
            end
            u = (1./nm(Nx)).*Rh_outlet^(2/3)*out_slope^0.5; % Normal velocity at the outlet
            if wse_dif > 0
                flow_dir = 1; % Flowing towards the outlet
            else
                flow_dir = -1; % Flowing to inside of the channel
            end
        end        
else
    error('HydroHP not coded for Darcy-Weisbach. Wait for the new version or change the method for Manning.')
end

% Discharge
q2(Nx,1) = q1(Nx,1)*u*flow_dir; % Area x Velocity

%%% State Space Format %%%
% dq/dt + dF/dx = S, we solve for A(x,t) and Q(x,t)
% q = [A Q]' = [q1 q2]'
% F [Q (Qv + gAybar]' = [q2 (q2^2)/q1 + g.q1.ybar]' = [f1 f2]'
% where ybar is the centroid depth from the top
% S = [0 gA(I0 - If)]'

% ybar = y - int(A(y)) / A(y) from y = 0 to y = y0
if flag_section ~= 4
    ybar = ybar_function(D,Z1,Z2,a,b,y0);
else
    % [y_irr, A_irr, P_irr, Rh_irr, y_bar_irr, n_med_irr, Beta_irr, u_irr, B_irr, Q_irr];
    % [   1,    2,     3,      4,    5,            6,          7,     8,      9,
    % ybar = y - ybar*
    %     ybar(:,1) =  Vlookup_leq(irr_table,col1,Q0*fv,1) - Vlookup_leq(irr_table,col1,Q0*fv,5);
    if (flag_hydrograph == 1 || flag_nash == 1) && flag_stage_hydrograph ~= 1
        col1 = 10; % Discharge
    else
        col1 = 1; % Depth
    end
    ybar(:,1) = Vlookup_g(irr_table,col1,s_v*fv,5);
end
f2(:,1) = q2(:,1).*abs(q2(:,1))./q1(:,1) + g*q1(:,1).*ybar(:,1);
f2(isnan(f2)) = 0; % Attention Here

% Friction S = [J1 J2]' with J1 = 0 and J2 calculated as follows:
if flag_friction == 1
    J2(:,1) = g*q1(:,1).*(I0(:) - q2(:,1).*abs(q2(:,1)).*nm(:).^2./(q1(:,1).^2.*Rh0.^(4/3))); % Manning
else
    J2(:,1) = g*q1(:,1).*(I0(:) - f*q2(:,1).*abs(q2(:,1))./((q1(:,1).^2).*8*g.*Rh0)); 
end

J2(isnan(J2)) = 0; % Attention Here

% Froude Number
if flag_section ~= 4
    Fr(:,1)=abs(q2(:,1)./q1(:,1))./((g*A_function(D,Z1,Z2,a,b,y0)./B_function(D,Z1,Z2,a,b,y0)).^0.5);% Froude Number
else
    % [y_irr, A_irr, P_irr, Rh_irr, y_bar_irr, n_med_irr, Beta_irr, u_irr, B_irr, Q_irr];
    % [   1,    2,     3,      4,    5,            6,          7,     8,      9,
    A_f_irr = Vlookup_g(irr_table,col1,s_v*fv,2)*ones(length(q1(:,1)),1);
    B_f_irr = Vlookup_g(irr_table,col1,s_v*fv,9)*ones(length(q1(:,1)),1);
    Fr(:,1)=abs(q2(:,1)./q1(:,1))./((g*A_f_irr./B_f_irr).^0.5);% Froude Number
end
% Courant Number
% Cn = c / (dx / dt), where c = v + sqrt(g.Hm), where Hm = A / B
if flag_section ~= 4
    Hm = A_function(D,Z1,Z2,a,b,y0)./B_function(D,Z1,Z2,a,b,y0);
    Cn(:,1)=(abs(q2(:,1)./q1(:,1))+(g*Hm).^0.5)/(dx/dt);% Courant Number
else
    Hm = A_f_irr./B_f_irr;
    Cn(:,1) = (abs(q2(:,1)./q1(:,1))+(g*Hm).^0.5)/(dx/dt);
end

% Depth in terms of Area function
% let c be the area in terms of Z1,Z2,b, and y, such that A(y) = c
% we want to solve y for A(y) = c

syms c_
if flag_section ~= 4
    fun_solve = (A - c_); % with c = area, we solve for y.
    options = optimoptions ('fsolve', 'Display', 'none','FunctionTolerance',1e-2,'MaxFunctionEvaluations',Nx*10);
end
if flag_section == 1
    % We have an analytical solution for this case
    z = solve(fun_solve,y_); % solving for y_ = y and c = A(y)
    h_function = matlabFunction(z); % h(A) = z;
else
    % Non-linear set of equations for circular pipe, we need to use fsolve
end
if flag_section ~= 4
    fun_solve = matlabFunction(fun_solve); % Transforming into an equation
end
%% 6.0 -  Main Loop %%
n = 1; % initializing counter
x_i = 2:(Nx-1); % vector for interior sections varying from 2 to (Nx - 1)
tic % starts measuring time
% Interpolation Variables
if flag_section == 4
    var_inlet = zeros((length(irr_table(1,:))-1),1,2); var_outlet = var_inlet;
    alfa_var_inlet = zeros((length(irr_table(1,:))-1),1,1); alfa_var_outlet = alfa_var_inlet;
    var_middle = zeros((length(irr_table(1,:))-1),length(x_i),2);
    alfa_var_middle = zeros((length(irr_table(1,:))-1),length(x_i),1);
end

% Initialization of some variables
time_end_min = (Nt)*dt;
time = 0;
time_previous = 0;
time_step = dt; % sec
t_store_prev = 0;

while time <= (time_end_min) % Main loop
    try
    n = n + 1; % Time-step index        
    time = time + time_step; % Seconds
    time_save_model(n) = time; % Seconds

    % Model Status
    percentage_timestepsec_maxCourant_maxh = [time/(tt)*100, time_step, max(max(Cn)), max(max(y))]

    % Agregating Inflows to the New Time-step
    if flag_hydrograph == 1 || flag_nash == 1
        z1 = find(tint > time_previous,1,'first'); % begin of the time-step
        z2 = find(tint <= time,1,'last'); % end of the time-step
        if isempty(z1)
            z1 = 1;
        end
        if isempty(z2) || z2 < z1
            z2 = z1;
        end
        if time_step >= dt
            Q0 = mean(Qe1int(z1:z2));
        else
            Q0 = Qe1int(z1);
        end
    end
    if time > 4.08*10^3
        ttt = 1;
    end
    % Agregating Stages to the New Time-step
    if flag_stage_hydrograph == 1
        z1 = find(tint > time_previous,1,'first'); % begin of the time-step
        z2 = find(tint <= time,1,'last'); % end of the time-step
        if isempty(z1)
            z1 = 1;
        end
        if isempty(z2) || z2 < z1
            z2 = z1;
        end
        if time_step >= dt
            h0 = mean(he1int(z1:z2));
        else
            h0 = he1int(z1);
        end
    end

    % Stop Program if Complex Number Occurs
    if imag(max(Cn(:,2-1))) > 0 || imag(max(q2(:,2-1)))
        error('Complex number possibly due to changing the regime from free flow to pressurized flow.')
    end
    %%%%% -  Boundary Conditions -  %%%%%
    %% Channel's begin (INLET)
    if flag_stage_hydrograph == 1
        %         h0 = he1int(n,1); % Water depth at x = 0 at time = time
        if flag_section == 4
            if h0 > max(irr_table(:,1))
                error('The maximum water depth is larger than the channel height.')
            end
            q1(1,2) = Vlookup_g(irr_table,1,h0,2); % Smaller values
        else
            q1(1,2) = A_function(D,Z1(Nx),Z2(Nx),a,b(Nx),h0);
        end
    else
        q1(1,2) = q1(2,1); % Area at section 1 is equals area of section 2 from previous time-step
    end
    if flag_hydrograph == 1 || flag_nash == 1
        %         q2(1,2) = Qe1int(n,1); % Flow at section 1 is the inflow hydrograph
        q2(1,2) = Q0; % Flow at section 1 is the inflow hydrograph
    else
        q2(1,2) = q2(2,1); % Flow at section 1 equals flow at section 2 from previous time-step
    end

    if flag_hydrograph == 0 && flag_nash == 0 && flag_stage_hydrograph == 0 && flag_outlet == 0
        q2(1,2) = q2(2,1); % Flow at section 1 equals flow at section 2 from previous time-step
%         q2(1,1) = q2(2,1);
    end

    % Interpolating All Values from I_rr_table using q1 as basis
    % Explanation: area is given in m2. P, Rh, and other variables are
    % in m. So we have a quadratically similar triangle relationship
    if flag_section == 4
        for mm = 1:(length(irr_table(1,:))-1)
            interp_base = q1(1,2); % Value that will be used for interpolation (area)
            if interp_base <= min_area % Col with area = 0
                area_smaller = 0; % Smaller values
            else
                area_smaller = Vlookup_l(irr_table,2,interp_base,2); % Smaller values
            end
            area_larger = Vlookup_g(irr_table,2,interp_base,2); % Larger values
            col1 = 2; % Interpolating from area values
            if interp_base <= min_area % Col with area = 0
                var_inlet(mm,1,1) = irr_table(2,mm); % Smaller values
            else
                var_inlet(mm,1,1) = Vlookup_l(irr_table,col1,interp_base,mm); % Smaller values
            end            
            var_inlet(mm,1,2) = Vlookup_g(irr_table,col1,interp_base,mm); % Larger values
            alfa_var_inlet(mm,1) = sqrt((interp_base - area_smaller)/(area_larger - area_smaller));
        end
    end

    if flag_section == 1 % Trapezoid or Rectangular
        if Z1(1) > 0 || Z2(1) > 0 % Trapezoidal channel
            y(1,2) = max(h_function(D,Z1(1),Z2(1),a,b(1),q1(1,2)')); % water depth in terms of area q1
            % In this previous function, we solve h = y in terms of A = q1 = c
        else
            y(1,2) = q1(1,2)/b; % water depth in terms of area q1 for rectangular channels
        end
    elseif flag_section > 1 % circular or paraboloid or irregular
        y0_guess = y(1,2-1);
        c = q1(1,2)*fv;  % WEIRDO. I HAVE TO CHECK IT OUT ... ISNT IT (2-1)?
        if flag_section ~= 4
            fun = @(y_)fun_solve(D,Z1(1),Z2(2),a,b(1),c,y_);
            y(1,2) = fsolve(fun,y0_guess,options); % non-linear solver
        else % Irregular section
            % [y_irr, A_irr, P_irr, Rh_irr, y_bar_irr, n_med_irr, Beta_irr, u_irr, B_irr, Q_irr];
            % [   1,    2,     3,      4,    5,            6,          7,     8,      9,
            col1 = 2; % Col with A
            col_var = 1;
            % Var* = Var(-) + alfa*(Var(+) - Var(-))
            y(1,2) = var_inlet(col_var,1,1) + alfa_var_inlet(col_var,1)*(var_inlet(col_var,1,2) - var_inlet(col_var,1,1));
            %             y(1,2) = Vlookup_leq(irr_table,col1,c,1);
        end
    end
    % ybar
    if flag_section ~= 4
        ybar(1,2) = ybar_function(D,Z1(1),Z2(1),a,b(1),y(1,2));
        % f1 and f2
        f1(1,2) = q2(1,2);
        f2(1,2) = q2(1,2).*abs(q2(1,2))./q1(1,2) + g*q1(1,2).*ybar(1,2);
        % Hydraulic Radius
        Rh_inlet = Rh_function(D,Z1(1),Z2(1),a,b(1),y(1,2));
        % Friction
        if flag_friction == 1
            J2(1,2) = g*q1(1,2).*(I0(1) - q2(1,2).*abs(q2(1,2)).*nm(1).^2./(q1(1,2).^2*Rh_inlet.^(4/3))); % Manning
        else
            J2(1,2) = (I0(1) - f*q2(1,2).*abs(q2(1,2))./((q1(1,2).^2)*8*g.*Rh_inlet));
        end
        % Froude
        Fr(1,2)=abs(q2(1,2)./q1(1,2))./((g*A_function(D,Z1(1),Z2(1),a,b(1),y(1,2))./B_function(D,Z1(1),Z2(1),a,b(1),y(1,2)))^0.5);% Froude Number
        % Courant
        Hm = A_function(D,Z1(1),Z2(1),a,b(1),y(1,2))./B_function(D,Z1(1),Z2(1),a,b(1),y(1,2));
        Cn(1,2)=(abs(q2(1,2)./q1(1,2))+(g*Hm).^0.5)/(dx/time_step);% Courant Number
        if isinf(Cn(1,2)) || isnan(Cn(1,2))
            Cn(1,2) = 0;
        end
    else
        % [y_irr, A_irr, P_irr, Rh_irr, y_bar_irr, n_med_irr, Beta_irr, u_irr, B_irr, Q_irr];
        % [   1,    2,     3,      4,    5,            6,          7,  8,      9,  10]
        col1 = 2; % Col with A
        col_var = 5;
        % Var* = Var(-) + alfa*(Var(+) - Var(-))
        ybar(1,2) = var_inlet(col_var,1,1) + alfa_var_inlet(col_var,1)*(var_inlet(col_var,1,2) - var_inlet(col_var,1,1));
        % f1 and f2
        f1(1,2) = q2(1,2);
        f2(1,2) = q2(1,2).*abs(q2(1,2))./q1(1,2) + g*q1(1,2).*ybar(1,2);
        % Hydraulic Radius
        col_var = 4;
        % Var* = Var(-) + alfa*(Var(+) - Var(-))
        Rh_inlet = var_inlet(col_var,1,1) + alfa_var_inlet(col_var,1)*(var_inlet(col_var,1,2) - var_inlet(col_var,1,1));
        % Friction
        if flag_friction == 1
            col_var = 6;
            % Var* = Var(-) + alfa*(Var(+) - Var(-))
            nm(1) = var_inlet(col_var,1,1) + alfa_var_inlet(col_var,1)*(var_inlet(col_var,1,2) - var_inlet(col_var,1,1));
            if isnan(nm(1,1))
                nm = irr_table(2,6)*ones(length(q1(:,1)),1);
            end
            J2(1,2) = g*c.*(I0(1) - q2(1,2).*abs(q2(1,2)).*nm(1).^2./(c.^2*Rh_inlet.^(4/3))); % Manning
        else
            J2(1,2) = (I0(1) - f*q2(1,2).*abs(q2(1,2))./((q1(1,2).^2)*8*g.*Rh_inlet));
        end
        % Froude
        % Var* = Var(-) + alfa*(Var(+) - Var(-))
        A_f_irr = q1(1,2);
        col_var = 9;
        B_f_irr = var_inlet(col_var,1,1) + alfa_var_inlet(col_var,1)*(var_inlet(col_var,1,2) - var_inlet(col_var,1,1));
        %         B_f_irr = Vlookup_leq(irr_table,col1,c,9);
        Fr(1,2) = abs(q2(1,2)./q1(1,2))./((g*A_f_irr./B_f_irr)^0.5);% Froude Number
        % Courant
        Hm = A_f_irr./B_f_irr;
        Cn(1,2) = (abs(q2(1,2)./q1(1,2))+(g*Hm).^0.5)/(dx/time_step);% Courant Number
    end

    %% Right side of the channel (outlet)
    if flag_outlet == 1 % Normal Depth
        q1(Nx,2) = q1(Nx-1,2-1); % Boundary Condition (same area)
        % Interpolating All Values from I_rr_table using q1 as basis
        % Explanation: area is given in m2. P, Rh, and other variables are
        % in m. So we have a quadratically similar triangle relationship
        if flag_section == 4
            for mm = 1:(length(irr_table(1,:))-1)
                interp_base = q1(Nx,2); % Value that will be used for interpolation (area)
                if interp_base <= min_area % Area
                    area_smaller = 0;
                else
                    area_smaller = Vlookup_l(irr_table,2,interp_base,2); % Smaller values
                end
                area_larger = Vlookup_g(irr_table,2,interp_base,2); % Larger values
                col1 = 2; % Interpolating from area values
                if interp_base <= min_area % Area
                    var_outlet(mm,1,1) = irr_table(2,mm); % Smaller values
                else
                    var_outlet(mm,1,1) = Vlookup_l(irr_table,col1,interp_base,mm); % Smaller values
                end
                var_outlet(mm,1,2) = Vlookup_g(irr_table,col1,interp_base,mm); % Larger values
                alfa_var_outlet(mm,1) = sqrt((interp_base - area_smaller)/(area_larger - area_smaller));
            end
        end
        if flag_section == 1
            if Z1(Nx) > 0 || Z2(Nx) > 0
                y(Nx,2) = max(h_function(D,Z1(Nx),Z2(Nx),a,b(Nx),q1(Nx,2)')); % water depth in terms of area q1
            else
                y(Nx,2) = q1(Nx,2)/b; % water depth in terms of area q1 for rectangular channels
            end
        elseif flag_section >= 2 % circular or paraboloid or irregular
            % If we do not have a stage-hydrograph boundary condition
            y0_guess = y(Nx,2-1);
            if flag_section ~= 4
                fun = @(y_)fun_solve(D,Z1(Nx),Z2(Nx),a,b(Nx),c,y_);
                y(Nx,2) = fsolve(fun,y0_guess,options); % non-linear solver
            else
                % [y_table, A, P, Rh, y_bar, n_med, Beta, v, B, Q]
                % [   1,    2, 3, 4,    5,    6,      7,  8,  9, 10]
                col_var = 1;
                % Var* = Var(-) + alfa*(Var(+) - Var(-))
                y(Nx,2) = var_outlet(col_var,1,1) + alfa_var_outlet(col_var,1)*(var_outlet(col_var,1,2) - var_outlet(col_var,1,1));
            end
        end
    else
        % Stage Hydrograph Boundary Condition. We are modeling a tidal
        % outlet condition
        time_wave = time; % time in seconds
        y(Nx,2) = h_wave_function(time_wave);
        if flag_section ~= 4
            q1(Nx,2) = A_function(D,Z1(Nx),Z2(Nx),a,b(Nx),y(Nx,2));
        else
            % We search Everything Using the Depth instead of the Flow
            col1 = 1; % Searching with the Col of Flow


            area_smaller = Vlookup_l(irr_table,col1,y(Nx,2)*fv,2);
            area_greater = Vlookup_g(irr_table,col1,y(Nx,2)*fv,2);
            y_smaller = Vlookup_l(irr_table,col1,y(Nx,2)*fv,1);
            y_greater = Vlookup_g(irr_table,col1,y(Nx,2)*fv,1);

            delta_y = y(Nx,2) - Vlookup_l(irr_table,col1,y(Nx,2)*fv,1);
            q1(Nx,2) = area_smaller + (area_greater - area_smaller)*(delta_y/(y_greater - y_smaller))^2;
        end
        %         q1(Nx,2) = q1(Nx-1,2-1)
    end
    % Hydraulic Radius
    if flag_section ~= 4
        Rh_outlet = Rh_function(D,Z1(Nx),Z2(Nx),a,b(Nx),y(Nx,2));
    else
        for mm = 1:(length(irr_table(1,:))-1)
            interp_base = q1(Nx,2); % Value that will be used for interpolation (area)
            if interp_base <= min_area
                area_smaller = 0; % Smaller values
            else
                area_smaller = Vlookup_l(irr_table,2,interp_base,2); % Smaller values
            end
            area_larger = Vlookup_g(irr_table,2,interp_base,2); % Larger values
            col1 = 2; % Interpolating from area values
            if interp_base <= min_area
                var_outlet(mm,1,1) = irr_table(2,mm);
            else
                var_outlet(mm,1,1) = Vlookup_l(irr_table,col1,interp_base,mm); % Smaller values
            end
            var_outlet(mm,1,2) = Vlookup_g(irr_table,col1,interp_base,mm); % Larger values
            alfa_var_outlet(mm,1) = sqrt((interp_base - area_smaller)/(area_larger - area_smaller));
        end
        % [y_table, A, P, Rh, y_bar, n_med, Beta, v, B, Q]
        % [   1,    2, 3, 4,    5,    6,      7,  8,  9, 10]
        col_var = 4; % Calculating Hydraulic Radius
        % Var* = Var(-) + alfa*(Var(+) - Var(-))
        Rh_outlet = var_outlet(col_var,1,1) + alfa_var_outlet(col_var,1)*(var_outlet(col_var,1,2) - var_outlet(col_var,1,1)); % Interpolated Hydraulic Radius
        % Var* = Var(-) + alfa*(Var(+) - Var(-))
        col_var = 6;
        nm(end,1) = var_outlet(col_var,1,1) + alfa_var_outlet(col_var,1)*(var_outlet(col_var,1,2) - var_outlet(col_var,1,1));
    end
    if flag_friction == 1
        if flag_outlet == 1
            u = (1./nm(Nx)).*Rh_outlet^(2/3)*I0(Nx)^0.5; % Normal depth at the outlet
            flow_dir = 1;
        else
            wse_dif = y(Nx-1,2-1) + inv_el(Nx-1) - y(Nx,2) - inv_el(Nx); % Difference in wse
            out_slope = abs(wse_dif)/dx; % Friction slope at the outlet as a diffusive model
            if wse_dif < 0 
                ttt = 1;
            end
              
            u = (1./nm(Nx)).*Rh_outlet^(2/3)*out_slope^0.5; % Normal velocity at the outlet
            if wse_dif > 0
                flow_dir = 1; % Flowing towards the outlet
            else
                flow_dir = -1; % Flowing to inside of the channel
            end
        end
    else
        u = sqrt(8*g*Rh_outlet*I0(Nx)/f); % outlet velocity
    end
    % Outlet Flow
    q2(Nx,2) = q1(Nx,2)*u*flow_dir; % Area x Velocity
    if isnan(q2(Nx,2))
        ttt = 1;
    end
    % Outlet Flow Under No Inflow Hydrograph & Not Enough WSE_dif
    if flag_stage_hydrograph ~= 1 && flag_nash ~= 1 && flag_hydrograph ~= 1
        % Only Outlet Tidal B.C.
        if wse_dif > 0 && y(Nx-1,2-1) <= fv*1e-3
            q2(Nx,2) = q1(Nx,2)*dx/(time_step); % Making sure all available depth becomes outflow in the outlet
        end
    end

    % ybar
    if flag_section ~= 4
        ybar(Nx,2) = ybar_function(D,Z1(Nx),Z2(Nx),a,b(Nx),y(Nx,2));
    else
        % [y_irr, A_irr, P_irr, Rh_irr, y_bar_irr, n_med_irr, Beta_irr, u_irr, B_irr, Q_irr];
        % [   1,    2,     3,      4,    5,            6,          7,  8,      9,  10]
        % ybar = y - ybar*
        col1 = 2; % A
        %         ybar(Nx,2) =  Vlookup_leq(irr_table,col1,c,1) - Vlookup_leq(irr_table,col1,c,5);
%         if q1(Nx,2) == 0
%             ybar(Nx,2) = 0;
%         else
%             ybar(Nx,2) = Vlookup_l(irr_table,col1,q1(Nx,2),5);
%         end
        col_var = 5;
        % Var* = Var(-) + alfa*(Var(+) - Var(-))
        ybar(Nx,2) = var_outlet(col_var,1,1) + alfa_var_outlet(col_var,1)*(var_outlet(col_var,1,2) - var_outlet(col_var,1,1));
    end
    % f1 and f2
    f1(Nx,2) = q2(Nx,2); % f1 - Flow
    zzz = q2(Nx,2).*abs(q2(Nx,2))./q1(Nx,2) + g*q1(Nx,2).*ybar(Nx,2); % f2 = (Qv + gAy_bar)
    zzz(isnan(zzz)) = 0;
    f2(Nx,2) = zzz; % f2 = (Qv + gAy_bar)

    % J2
    % Friction
    if flag_friction == 1
        J2(Nx,2) = g*q1(Nx,2).*(I0(Nx) - q2(Nx,2).*abs(q2(Nx,2)).*nm(Nx)^2./(q1(Nx,2).^2*Rh_outlet.^(4/3))); % Manning --> gA*(I0 - If), If = n^2*Q*abs*Q)/(Rh^(4/3)*A^2)
    else
        J2(Nx,2) = g*q1(Nx,2).*(I0(Nx) - f*q2(:,2).*abs(q2(Nx,2))./((q1(Nx,2).^2)*8*g*Rh_outlet));
    end
    J2(isnan(J2)) = 0; % Attention Here
    % Froude
    if flag_section ~= 4
        Fr(Nx,2)=abs(q2(Nx,2)./q1(Nx,2))./((g*A_function(D,Z1(Nx),Z2(Nx),a,b(Nx),y(Nx,2))./B_function(D,Z1(Nx),Z2(Nx),a,b(Nx),y(Nx,2)))^0.5);% Froude Number
        % Courant
        Hm = A_function(D,Z1(Nx),Z2(Nx),a,b(Nx),y(Nx,2))./B_function(D,Z1(Nx),Z2(Nx),a,b(Nx),y(Nx,2));
        Cn(Nx,2)=(abs(q2(Nx,2)./q1(Nx,2))+(g*Hm).^0.5)/(dx/time_step);% Courant Number
        if isnan(Cn(Nx,2)) || isinf(Cn(Nx,2))
            Cn(Nx,2) = 0;
        end
    else
        % Froude
        % [y_irr, A_irr, P_irr, Rh_irr, y_bar_irr, n_med_irr, Beta_irr, u_irr, B_irr, Q_irr];
        % [   1,    2,     3,      4,    5,            6,          7,  8,      9,  10]
        col1 = 2; % Col with A
        A_f_irr = c;
        col_var = 9;
        % Var* = Var(-) + alfa*(Var(+) - Var(-))
        B_f_irr = var_outlet(col_var,1,1) + alfa_var_outlet(col_var,1)*(var_outlet(col_var,1,2) - var_outlet(col_var,1,1));
        %         B_f_irr = Vlookup_l(irr_table,col1,c,9);
        Fr(Nx,2) = abs(q2(Nx,2)./q1(Nx,2))./((g*A_f_irr./B_f_irr)^0.5);% Froude Number
        % Courant
        Hm = A_f_irr./B_f_irr;
        if y(Nx,2) <= min_depth
            Cn(Nx,2) = 0;
        else
            Cn(Nx,2)=(abs(q2(Nx,2)./q1(Nx,2))+(g*Hm).^0.5)/(dx/time_step);% Courant Number
            if isnan(Cn(Nx,2)) || isinf(Cn(Nx,2))
                Cn(Nx,2) = 0;
            end   
        end
    end

    %% Main Loop for Non-Boundary Cells from 2 to (Nx - 1)
    % vectorized calculations
    q1_back = q1(1:(Nx-2),(2-1));
    q1_forward = q1(3:(Nx),(2-1));
    q2_back = q2(1:(Nx-2),(2-1));
    q2_forward = q2(3:(Nx),(2-1));
    f1_back = f1(1:(Nx-2),(2-1));
    f1_forward = f1(3:(Nx),(2-1));
    f2_back = f2(1:(Nx-2),(2-1));
    f2_forward = f2(3:(Nx),(2-1));
    J2_back = J2(1:(Nx-2),(2-1));
    J2_forward = J2(3:(Nx),(2-1));

    % Lax-Friedrichs Method
    % Given a hyperbolic partial derivative system of equations described
    % by:
    % pq/pt + pF/px - S = 0, where p is the partial derivative, one can
    % solve this equation by performing a forward discretization for q and a
    % central discretization for F. Moreover, S = (Sback + Sforward)/2
    % Expliciting the system of equations for q, it follows that:

    q1(x_i,2) = 0.5.*(q1_forward + q1_back) - 0.5*time_step/dx*(f1_forward - f1_back); %% attention here in f1forward
    q2(x_i,2) = 0.5*(q2_forward + q2_back) - 0.5*time_step/dx*(f2_forward - f2_back) + 0.5*time_step*(J2_back + J2_forward);
    
    if q1(Nx-1,2) > 0.0
        ttt = 1;
    end
    % There is no such thing as a negative water depth, so we apply a
    % constraint
%     if min(q1(x_i,2)) < 0
%         zzz = q1(x_i,2); zzz(zzz<0) = 0; q1(x_i,2) = zzz;
%         ttt = 1;
%     end
    % Interpolating All Values from I_rr_table using q1 as basis
    if flag_section == 4
        for mm = 1:(length(irr_table(1,:))-1)
            for hh = 1:length(x_i)
                interp_base = q1(hh+1,2); % Value that will be used for interpolation (area)
                if interp_base <= min_area
                    area_smaller = 0;
                else
                    area_smaller = Vlookup_l(irr_table,2,interp_base,2); % Smaller values   
                end
                area_larger = Vlookup_g(irr_table,2,interp_base,2); % Larger values
                col1 = 2; % Interpolating from area values
                if interp_base <= min_area
                    var_middle(mm,hh,1) = irr_table(2,mm); % Smaller values
                else
                    var_middle(mm,hh,1) = Vlookup_l(irr_table,col1,interp_base,mm); % Smaller values
                end
                var_middle(mm,hh,2) = Vlookup_g(irr_table,col1,interp_base,mm); % Larger values
                alfa_var_middle(mm,hh,1) = sqrt((interp_base - area_smaller)/(area_larger - area_smaller));
            end
        end
    end

    if flag_section == 1
        if Z1(2) >0 || Z2(2) >0
            y(x_i,2) = max(h_function(D,Z1(x_i),Z2(x_i),a,b(x_i),q1(x_i,2)')); % water depth in terms of area q1
        else
            y(x_i,2)=q1(x_i,2)/b;
        end
    elseif flag_section > 1
        y0_guess = y(x_i,2-1);
        c = q1(x_i,2)*fv; % It has to be a line vector (area)
        if flag_section ~= 4
            fun = @(y_)fun_solve(D,Z1(x_i),Z2(x_i),a,b(x_i),c,y_);
            y(x_i,2) = fsolve(fun,y0_guess,options); % non-linear solver
        else
            % [y_irr, A_irr, P_irr, Rh_irr, y_bar_irr, n_med_irr, Beta_irr, u_irr, B_irr, Q_irr];
            % [   1,    2,     3,      4,    5,            6,          7,  8,      9,  10]
            col1 = 2; % Col with A
            for i = 1:length(x_i)
                cc = c(i); % be careful here
                col_var = 1;
                % Var* = Var(-) + alfa*(Var(+) - Var(-))
                y(i+1,2) = var_middle(col_var,i,1) + alfa_var_middle(col_var,i)*(var_middle(col_var,i,2) - var_middle(col_var,i,1));
            end
        end
    end
    % Hydraulic Radius
    if flag_section ~= 4
        Rh_middle = Rh_function(D,Z1(x_i),Z2(x_i),a,b(x_i),y(x_i,2));
        % ybar
        ybar(x_i,2) = ybar_function(D,Z1(x_i),Z2(x_i),a,b(x_i),y(x_i,2));
        % f1 and f2
        f1(x_i,2) = q2(x_i,2);
        f2(x_i,2) = q2(x_i,2).*abs(q2(x_i,2))./q1(x_i,2) + g*q1(x_i,2).*ybar(x_i,2);
        % Froude
        Hm = A_function(D,Z1(x_i),Z2(x_i),a,b(x_i),y(x_i,2))./B_function(D,Z1(x_i),Z2(x_i),a,b(x_i),y(x_i,2));
        Fr(x_i,2)=abs(q2(x_i,2)./q1(x_i,2))./((g*Hm).^0.5);% Froude Number
        % Courant
        Cn(x_i,2)=(abs(q2(x_i,2)./q1(x_i,2))+(g*Hm).^0.5)/(dx/time_step);% Courant Number
        % Friction
        if flag_friction == 1
            J2(x_i,2) = g*q1(x_i,2).*(I0(x_i) - q2(x_i,2).*abs(q2(x_i,2).*nm(x_i).^2./(q1(x_i,2).^2.*Rh_middle.^(4/3))));
        else
            J2(x_i,2) = g*q1(x_i,2).*(I0(x_i) - f*q2(x_i,2).*abs(q2(x_i,2))./((q1(x_i,2).^2)*8*g*Rh_midle));
        end
        % Stability Check
        if max(Cn(:,2)) > 1
            error('Please, decrease the time-step')
        end
    else
        for jj = 1:length(x_i)
            cc = c(jj); % Area
            % [y_irr, A_irr, P_irr, Rh_irr, y_bar_irr, n_med_irr, Beta_irr, u_irr, B_irr, Q_irr];
            % [   1,    2,     3,      4,    5,            6,          7,  8,      9,  10]
            col_var = 4;
            % Var* = Var(-) + alfa*(Var(+) - Var(-))
            Rh_middle(jj,1) = var_middle(col_var,jj,1) + alfa_var_middle(col_var,jj)*(var_middle(col_var,jj,2) - var_middle(col_var,jj,1));
            col_var = 5;
            ybar(jj+1,2) = var_middle(col_var,jj,1) + alfa_var_middle(col_var,jj)*(var_middle(col_var,jj,2) - var_middle(col_var,jj,1));
            col_var = 6;
            nm(jj+1,1) = var_middle(col_var,jj,1) + alfa_var_middle(col_var,jj)*(var_middle(col_var,jj,2) - var_middle(col_var,jj,1));
            % f1 and f2
            f1(jj+1,2) = q2(jj+1,2);
            f2(jj+1,2) = q2(jj+1,2).*abs(q2(jj+1,2))./q1(jj+1,2) + g*q1(jj+1,2).*ybar(jj+1,2);
            % Froude
            A_f_irr = q1(jj+1,2);
            col_var = 9;
            B_f_irr = var_middle(col_var,jj,1) + alfa_var_middle(col_var,jj)*(var_middle(col_var,jj,2) - var_middle(col_var,jj,1));
            Hm = A_f_irr./B_f_irr;
            Fr(jj+1,2) = abs(q2(jj+1,2)./q1(jj+1,2))./((g*Hm).^0.5);% Froude Number
            % Courant
            if y(jj+1,2) > 0.005 % 0.5 cm
                Cn(jj+1,2) = (abs(q2(jj+1,2)./q1(jj+1,2))+(g*Hm).^0.5)/(dx/time_step);% Courant Number                
            else
                Cn(jj+1,2) = 0;
            end
            if isinf(Cn(jj+1,2))
                Cn(jj+1,2) = 0;
            end
            % Friction
            if flag_friction == 1
                J2(jj+1,2) = g*A_f_irr.*(I0(jj+1,1) - q2(jj+1,2).*abs(q2(jj+1,2).*nm(jj+1,1).^2./(A_f_irr.^2.*Rh_middle(jj,1)^(4/3))));
            else
                J2(jj+1,2) = g*q1(jj+1,2).*(I0(jj+1,2) - f*q2(jj+1,2).*abs(q2(jj+1,2))./((q1(jj+1,2).^2)*8*g*Rh_midle(jj,1)));
            end
            % Stability Check
            if Cn(jj+1,2) > 1 && y(jj+1,2) >= min_depth && q1(jj+1,2) >= min_area
                error('Please, decrease the time-step')
            end
        end
    end


    % Constraint at dry areas
    % -- the idea is that dry cells have no hydraulic properties
    if min(q1(:,2)) <=  min_area || min(y(:,2)) <=  min_depth
        idx1 = q1(:,2) <= min_area; idx2 = y(:,2) <= min_area; idx = idx1 + idx2; % Both
        idx = logical([zeros(size(idx,1),1), idx]);
        q1(idx) = 0; q2(idx) = 0; f1(idx) = 0; f2(idx) = 0; J2(idx) = 0; y(idx) = 0;
        Fr(idx) = 0; Cn(idx) = 0; 
    end
        % Adaptive Time-Step - Outlet not considered
        idx_courant = Cn(1:end-1,2) <= 0;
        zzz = Cn(1:end-1,2); zzz(idx_courant) = nan;
        dt_courant_1 = zzz/time_step; % Cn/time_step = dx/ (v + sqrt(Hm*g)), this the the time-step for Courant = 1
        time_step = min(alpha./dt_courant_1); % Calculated
        time_step = min(time_step,dtmax);
        time_step = max(time_step,dtmin);
        time_previous = time;    


    % Saving hydrographs and depths with user defined recording time-step
    if n == 1
        % Do nothing, it is already solved, we just have to save the data
        % for the next time-step
        t_store = 1;
        time_save(1,1) = time;
    else
        t_store = find(time_store <= time,1,'last'); % Time that is being recorded in min
        if t_store > t_store_prev
            time_save(t_store,1) = time;
            Flow_Area(t_store,:) = q1(:,2); % m2
            Discharge(t_store,:) = q2(:,2); % m3/s
            Depth(t_store,:) = y(:,2); % m  
            Velocity(t_store,:) = q2(:,2)./q1(:,1); % m/s
            Froude(t_store,:) = Fr(:,2);
            Courant(t_store,:) = Cn(:,2);
            t_store_prev = t_store;
        end
    end
    % Refreshing States
%     idx = y < 1e-3; q1(idx) = 0; q2(idx) 

    q1(:,1) = q1(:,2);
    q2(:,1) = q2(:,2);
    f1(:,1) = f1(:,2);
    f2(:,1) = f2(:,2);
    J2(:,1) = J2(:,2);
    y(:,1) = y(:,2);
    Cn(:,1) = Cn(:,2);
    Fr(:,1) = Fr(:,2);

    if time > 2*1000
        ttt = 1;
    end
    
    catch ME
        % If this condition is reached, we are reducing the time-step to
        % 50% and doing the calculations again
        idx = q1(:,2) <= min_area;
        vel = abs(q2(:,2)./q1(:,2)) + sqrt(g*y(:,2)); vel(idx) = 0;
        dtnew = min(alpha*dx./(vel));
        time = time - time_step; % Seconds
        time_step = dtnew; % Halving the time-step
        n = n - 1;        
    end
end
%% 7.0 - Post-Processing
water_depths = Depth;
%%% Post Processing Figures %%%
% Call function
warning('on');
post_processing
close all
computational_time
disp(['Thank you for using HydroHP. If you have any questions, please contact me at marcusnobrega.engcivil@gmail.com).'])
disp(['Also, please check your current matlab folder. The outputs are there.'])







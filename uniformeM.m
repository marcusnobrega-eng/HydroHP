%%%% Newton Raphson - Numerical Scheme

% Manning Equation
% Q = 1/n x A x (A / P) ^ (2/3) x (I0)^0.5
% Q x n x P^(2/3) = A ^(5/3) x I0^(0.5) = Fi

function y0 = uniformeM(nm,Q0,b,Z1,Z2,a,D,I0,P,A,y0_guess) 
syms D_ I0_ Q_ Z1_ Z2_ a_ b_ n_ y_
%% Fi
Fi = n_.*Q_.*P.^(2/3) - (I0_.^0.5).*(A).^(5/3);
Fi_function = matlabFunction(Fi); % Function describring Fi
%% dFi
dFi = diff(Fi,y_); % First order derivative of Fi
dFi_function = matlabFunction(dFi); % Function describing the derivative of Fi
%% Newton Raphson Method
%y0(k+1) = y0(k) - Fi / (dFi / dy)
k = 1;
y0N(1,1) = y0_guess; % initial guess of water level in meters
ep = 10^3;
y0 = zeros(length(I0),1);
for i=1:length(I0)
    while ep > 1e-4
        % y0(k+1) = y0(k) - Fi / (dFi / dy)
        % Calling Fi Function
        Fi = Fi_function(D,I0(i),Q0,Z1(i),Z2(i),a,b(i),nm(i),y0N(k,1)); % Please, make sure you are inputing the right paramters at this function in the right order
        if imag(Fi) ~= 0
            y0N(1,1) = y0N(1,1)/2; % half of initial guess
           % Calling Fi Function
            Fi = Fi_function(D,I0(i),Q0,Z1(i),Z2(i),a,b(i),nm(i),y0N(k,1)); % Please, make sure you are inputing the right paramters at this function in the right order
            if imag(Fi) ~= 0
                y0N(1,1) = y0N(1,1)/2; % half of initial guess
                % Calling Fi Function
                Fi = Fi_function(D,I0(i),Q0,Z1(i),Z2(i),a,b(i),nm(i),y0N(k,1)); % Please, make sure you are inputing the right paramters at this function in the right order
                if imag(Fi) ~= 0
                    y0N(1,1) = 0.01; % 1 cm minimum tolerance
                end
            end
        end
        dFi = dFi_function(D,I0(i),Q0,Z1(i),Z2(i),a,b(i),nm(i),y0N(k,1));
        y0N(k+1,1) = y0N(k) - Fi/dFi;
        ep = abs(y0N(k+1,1) - y0N(k,1)); % error
        k = k + 1; % increase counter
    end
    y0(i) = y0N(k-1,1); % final water depth in meters
    ep = 10^3; % Returning the error to the initial error
end
end
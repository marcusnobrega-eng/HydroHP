% Objective: Interpolate Irregular Trapezoid Cross-Sections
% Developer: Marcus Nobrega, Ph.D.
% Input_Data
% b1, Zl1, Zr1
% b2, Zl2, Zr2
% L = Distance within xsections
% x1 = Coordenate of the i section
% dx = Model node discretization

function [x_section,b_section,Zl_section,Zr_section,table] = generate_parameters_section(b1,Zl1,Zr1,b2,Zl2,Zr2,x1,dx,L)
   n_interval = floor((L)/dx);
   x_section = (1:1:n_interval)*dx - dx + x1;
   x_section = x_section';
   % Slopes
   alpha_b = (b2 - b1)/L;
   alpha_Zl = (Zl2 - Zl1)/L;
   alpha_Zr = (Zr2 - Zr1)/L;

   % Preallocation
   b_section = zeros(n_interval,1);
   Zl_section = zeros(n_interval,1);
   Zr_section = zeros(n_interval,1);

   for i = 1:n_interval
       % b
       b_section(i,1) = b1 + alpha_b*(i-1)*dx;
       % Zl
       Zl_section(i,1) = Zl1 + alpha_Zl*(i-1)*dx;
       % Zr
       Zr_section(i,1) = Zr1 + alpha_Zr*(i-1)*dx;       
   end

   % table
   table = [x_section, b_section, Zl_section, Zr_section];
end
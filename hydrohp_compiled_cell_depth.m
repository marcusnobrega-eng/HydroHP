function depth_m = hydrohp_compiled_cell_depth(cells, area_m2)
%HYDROHP_COMPILED_CELL_DEPTH Cell depth from area using compiled tables.
%
% cells is the 'cells' set of HYDROHP_COMPILE_GEOMETRY. The result equals
% interp1(table.area_m2, table.depth_m, area, 'pchip') for each cell, the
% value HYDROHP_DEPTH_FROM_AREA returns, with areas clamped to the table.

area_m2 = min(area_m2(:), cells.maximum_area_m2);
count = cells.count;
index = (1:count).';
low = ones(count, 1);
high = cells.rows;
while true
    open = high - low > 1;
    if ~any(open)
        break
    end
    middle = floor((low + high) / 2);
    values = cells.area_m2(index + (middle - 1) * count);
    go_right = open & values(:) <= area_m2;
    go_left = open & ~go_right;
    low(go_right) = middle(go_right);
    high(go_left) = middle(go_left);
end
base = index + (low - 1) * count;
local = area_m2 - cells.area_m2(base);
coefficients = cells.depth_of_area;
plane = count * size(coefficients, 2);
depth_m = ((coefficients(base) .* local + coefficients(base + plane)) .* local + ...
    coefficients(base + 2 * plane)) .* local + coefficients(base + 3 * plane);
end

function [output] = inbound(cellcoord)

% check if cell is inside domain
output = all([max(cellcoord)<[170,800],min(cellcoord)>[1,10]]);

end


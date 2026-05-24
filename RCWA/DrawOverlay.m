function DrawOverlay(sim, ~)
%ABSORPTION_SPECTRUM Summary of this function goes here
%   Detailed explanation goes here

Lx = sim.input.nmLx;
ztop = -(sim.setup.nmz(1) + sim.setup.nmdz(1)/2);
zbottom = -(sim.setup.nmz(end) - sim.setup.nmdz(end)/2);

line([-Lx/2, -Lx/2],[ztop, zbottom],'color',[0,0,0],'LineWidth',1.2); % Left Boundary
line([Lx/2, Lx/2],[ztop, zbottom],'color',[0,0,0],'LineWidth',1.2); % Right Boundary
line([-Lx/2, Lx/2],[ztop, ztop],'color',[0,0,0],'LineWidth',1.2); % Top Boundary
line([-Lx/2, Lx/2],[zbottom, zbottom],'color',[0,0,0],'LineWidth',1.2); % Bottom Boundary

zcurrent = zbottom;
for i = 1:length(sim.material.nmSec)
    zcurrent = zcurrent + sim.material.nmSec(end+1-i);
    line([-Lx/2, Lx/2],[zcurrent, zcurrent],'color',[0,0,0],'LineWidth',0.1); % Top Boundary
end

zcurrent = zbottom;
for i = 1:sim.setup.Nz
    
    nmdz = sim.setup.nmdz(end+1-i);
    for j = 1:size(sim.z.zetalist,1)
        zeta = sim.z.zetalist(j,end+1-i);
        line([zeta * Lx, zeta * Lx],[zcurrent-nmdz/2, zcurrent+nmdz/2],'color',[0,0,0],'LineWidth',2); % Left Grating
    end
    zcurrent = zcurrent + nmdz;
    % line([-zeta(j) * Lx/2, -zeta(j+1) * Lx/2],[zcurrent, zcurrent],'color',[0,0,0],'LineWidth',1); % Left Grating
    %line([zeta(j) * Lx/2, zeta(j+1) * Lx/2],[zcurrent, zcurrent],'color',[0,0,0],'LineWidth',1); % Left Grating
    % line([-zeta(j) * Lx/2, -zeta(j) * Lx/2],[zcurrent, zcurrent + nmdz(j) ],'color',[0,0,0],'LineWidth',1); % Left Grating
    % line([zeta(j) * Lx/2, zeta(j) * Lx/2],[zcurrent, zcurrent + nmdz(j) ],'color',[0,0,0],'LineWidth',1); % Right Grating
    % line([-zeta(end) * Lx/2, -zeta(end) * Lx/2],[zcurrent, zcurrent + nmdz(end) ],'color',[0,0,0],'LineWidth',1); % Left Grating
    % line([zeta(end) * Lx/2, zeta(end) * Lx/2],[zcurrent, zcurrent + nmdz(end) ],'color',[0,0,0],'LineWidth',1); % Right Grating
    
end



end


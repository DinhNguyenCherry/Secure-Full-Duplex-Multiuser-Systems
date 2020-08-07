function [ done ] = Plot_Layout( RadiusOfCell, InnerZoneRadius, positionDownlinkUsers, positionUplinkUsers, positionEavs )
%PLOT_LAYOUT Summary of this function goes here
%   Detailed explanation goes here

hold on

ang=0:0.01:2*pi; 
xp=RadiusOfCell*cos(ang);
yp=RadiusOfCell*sin(ang);
ax = plot(0+xp,0+yp);

xpi=InnerZoneRadius*cos(ang);
ypi=InnerZoneRadius*sin(ang);
plot(0+xpi,0+ypi);

sz = 100

scatter(0,0, 2*sz, 'k^ ')
scatter(positionUplinkUsers(:,1), positionUplinkUsers(:,2), sz, 'bs ')
scatter(positionDownlinkUsers(:,1), positionDownlinkUsers(:,2), sz, 'r* ')
scatter(positionEavs(:,1), positionEavs(:,2), sz, 'gx ')


set(gca,'XTick',[-100 :20: 100])
set(gca,'YTick',[-100 :20: 100])

xlim([-RadiusOfCell RadiusOfCell])
ylim([-RadiusOfCell RadiusOfCell])
axis('square')

done = 1;


end


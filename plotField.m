function [] = plotField(I, figNum)
%PLOTFIELD Plots the initial setup for an EKF run with robot position and 
%orientation and feature positions 
xr = I.x_init;
qr = I.q_init;
T_i2b = quat2dcm(qr');
map = I.map;
tags = I.ptTags;
figure(figNum);
hold on
for ii = 1:size(map,2)
    plot3(map(1,ii),map(2,ii),map(3,ii),'*r');
    text(map(1,ii)+.05,map(2,ii)+.05,map(3,ii)+.05,num2str(tags(ii)))
end  
plot3(xr(1),xr(2),xr(3),'*b')
axisParam = I.plotAxisParam;
xb = T_i2b'*[axisParam;0;0]; % active rotation b2i
yb = T_i2b'*[0;axisParam;0];
zb = T_i2b'*[0;0;axisParam];
quiver3(xr(1),xr(2),xr(3),xb(1),xb(2),xb(3),'r','LineWidth',2)
quiver3(xr(1),xr(2),xr(3),yb(1),yb(2),yb(3),'g','LineWidth',2)
quiver3(xr(1),xr(2),xr(3),zb(1),zb(2),zb(3),'b','LineWidth',2)
hold off
axis equal
view(45,45)
end


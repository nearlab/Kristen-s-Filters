function [] = plotTrajectory(tHist,statesHist,figNum,I)
%PLOTTRAJECTORY Plot a trajectory in figure figNum
lenT = length(tHist);
xHist = statesHist(:,1:3);
qHist = statesHist(:,7:10);
dt = floor(lenT/20);

% Plot trajectory
figure(figNum);
hold on
plot3(xHist(:,1),xHist(:,2),xHist(:,3));
axisParam = I.plotAxisParam;

% Plot 10 attitudes
for ii = 1:dt:lenT 
    xr = xHist(ii,:)';
    qr = qHist(ii,:)';
    T_i2b = quat2dcm(qr');
    xb = T_i2b'*[axisParam;0;0]; 
    yb = T_i2b'*[0;axisParam;0];
    zb = T_i2b'*[0;0;axisParam];
    quiver3(xr(1),xr(2),xr(3),xb(1),xb(2),xb(3),'r','LineWidth',2)
    quiver3(xr(1),xr(2),xr(3),yb(1),yb(2),yb(3),'g','LineWidth',2)
    quiver3(xr(1),xr(2),xr(3),zb(1),zb(2),zb(3),'b','LineWidth',2)
%     pause(0.1)
end
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
axis equal
end


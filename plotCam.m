function [] = plotCam(cm,mask,figNum)
%PLOTCAM Constructs the image taken by a camera.
%   cm is a 2xN vector of (x,y) feature positions in the image
%   mask removes the features not in view of the camera
figure(figNum); clf
plot(cm(1,:),cm(2,:),'*r')
axis([-2,2,-2,2])
pause(0.1)
end
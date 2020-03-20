function [] = plotCam(cm,pfTagsObserved,fov,figNum)
%PLOTCAM Constructs the image taken by a camera.
%   cm is a 2xN vector of (x,y) feature positions in the image
%   mask removes the features not in view of the camera
figure(figNum); clf
for ii = 1:size(cm,2)
    hold on; plot(cm(1,ii),cm(2,ii),'*r')
    text(cm(1,ii)+.05,cm(2,ii)+.05,num2str(pfTagsObserved(ii)))
end
rectangle('Position',[-fov,-fov,2*fov,2*fov])
axis([-2,2,-2,2])
pause(0.1)
end
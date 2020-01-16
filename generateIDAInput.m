function I = generateIDAInput()
%GENERATEHESCHINPUT Generates Hesch trajectory, initial states, and true
%states, and map
I.x_init = [0, 0, 0]';
I.v_init = [0.1, 0.1, 0]';
I.q_init = [0, 0, 0.7071, 0.7071]';
I.T = 10;
I.dt = 0.1;
I.tVec = 0:I.dt:I.T;
I.trueInput = @trueInput;
nCirclePts = 8;
mapCircle = [0.5+0.1*cos((1:nCirclePts).*(2*pi/nCirclePts)); 
             0.6*ones(1,nCirclePts);
             0.1*sin((1:nCirclePts).*(2*pi/nCirclePts))];
mapPolygon = [0.8, 0.6, 0.05;
              0.8, 0.6, -0.05;
              0.55, 0.6, -0.25;
              0.45, 0.6, -0.25;
              0.2, 0.6, 0.05;
              0.2, 0.6, -0.05;
              0.55, 0.6, 0.25;
              0.45, 0.6, 0.25;]';
I.map = [mapCircle, mapPolygon];
I.nPts = size(I.map,2);
I.ptTags = 1:size(I.map,2);
I.plotAxisParam = 0.1;
end

function [at_i,wt_b] = trueInput(t)
at_i = [-0.01;
        -0.01;
        0];
wt_b = [0;
        0;
        0];
end


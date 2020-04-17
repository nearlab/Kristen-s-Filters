function I = generateHeschInput()
%GENERATEHESCHINPUT Generates Hesch trajectory, initial states, and true
%states, and map
I.x_init = [5, 0, 0]';
I.v_init = [0, 2*pi*5/60, 0]';
I.q_init = [-0.5, -0.5, -0.5, -0.5]';
I.T = 120;
I.dt = 0.05;
I.tVec = 0:I.dt:I.T;
I.trueInput = @trueInput;
mapFreq = 2*pi/50;
mapLo = -0.25;
mapHi = 0.25;
mapHeight = 2;
[z,~] = meshgrid(linspace(mapLo,mapHi,mapHeight),mapFreq:mapFreq:2*pi);
I.map = [repmat(6*cos(mapFreq:mapFreq:2*pi),1,mapHeight);
         repmat(6*sin(mapFreq:mapFreq:2*pi),1,mapHeight);
         z(:)']; 
I.nPts = size(I.map,2);
I.ptTags = 1:size(I.map,2);
I.plotAxisParam = 1; % [m] length of axes in attitude plots
end

function [at_i,wt_b] = trueInput(t)
at_i = [-5*(2*pi/60)^2*cos((2*pi/60)*t);
        -5*(2*pi/60)^2*sin((2*pi/60)*t);
        0];
wt_b = [0;
        (2*pi/60);
        0];
end


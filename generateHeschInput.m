function I = generateHeschInput()
%GENERATEHESCHINPUT Generates Hesch trajectory, initial states, and true
%states, and map
I.x_init = [5, 0, 0]';
I.v_init = [0, 2*pi*5/60, 0]';
I.q_init = [-0.5, -0.5, -0.5, -0.5]';
I.T = 60;
I.dt = 0.25;
I.tVec = 0:I.dt:I.T;
I.trueInput = @trueInput;
mapFreq = 2*pi/36;
mapLo = -1;
mapHi = 1;
mapHeight = 6;
[z,~] = meshgrid(linspace(mapLo,mapHi,mapHeight),mapFreq:mapFreq:2*pi);
I.map = [repmat(6*cos(mapFreq:mapFreq:2*pi),1,mapHeight);
         repmat(6*sin(mapFreq:mapFreq:2*pi),1,6);
         z(:)']; 
I.nPts = size(I.map,2);
I.ptTags = 1:size(I.map,2);
I.plotAxisParam = 1; % [m] length of axes in attitude plots
end

function [at_i,wt_b] = trueInput(t)
at_i = [-5*(2*pi/60)^2*cos((2*pi/60)*t);
        -5*(2*pi/60)^2*sin((2*pi/60)*t);
        -1*(2*pi/10)^2*sin((2*pi/10)*t)];
wt_b = [zeros(1,length(t));
        zeros(1,length(t)); %(2*pi/60)*ones(1,length(t));
        zeros(1,length(t))];
end


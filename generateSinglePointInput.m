function I = generateSinglePointInput()
%GENERATEHESCHINPUT Generates Hesch trajectory, initial states, and true
%states, and map
I.x_init = [5, 0, 0]';
I.v_init = [0, 2*pi*5/60, 2*pi*1/10]';
I.q_init = [-0.5, -0.5, -0.5, -0.5]';
I.T = 60;
I.dt = 0.01;
I.tVec = 0:I.dt:I.T;
I.trueInput = @trueInput;
I.map = [6; 0; 0];
I.nPts = 1;
I.ptTags = 1;
I.plotAxisParam = 1; % [m] length of axes in attitude plots
end

function [at_i,wt_b,wt_i] = trueInput(t)
at_i = [-5*(2*pi/60)^2*cos((2*pi/60)*t);
        -5*(2*pi/60)^2*sin((2*pi/60)*t);
        -1*(2*pi/10)^2*sin((2*pi/10)*t)];
wt_b = [zeros(1,length(t));
        0*(2*pi/60)*ones(1,length(t));
        zeros(1,length(t))];
wt_i = [zeros(1,length(t));
        zeros(1,length(t));
        0*(2*pi/60)*ones(1,length(t));];
% wt_b = [0;
%         (2*pi/60);
%         0];
end


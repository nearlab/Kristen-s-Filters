%%%Inputs
close all;
I = generateSinglePointInput();
qf = quaternionFunctions();

%initial state
x_init = I.x_init;
v_init = I.v_init;  %inertial frame
q_init = I.q_init;
T = I.T; %must be an even multiple of dt
dt = I.dt;
tVec = I.tVec;
plotField(I,1);
tspan = [0,T];

states_init = [x_init; v_init; q_init];
[t_hist,states_hist_true] = ode45(@(t,states)trueKinematics(t,states,I,qf),tVec,states_init);
plotTrajectory(t_hist,states_hist_true,1,I)
%uncertainties
%x: +/- .01 m
%v: +/- .01 m/s
%q: +/- .1 rad
PSDa = (.7/60)^2; % [(m/s)^2/s] = 0.7 m/s/sqrt(hr)
PSDg = (0.15*pi/180/60)^2; % [rad^2/s] = 0.15 deg/sqrt(hr)
Px = diag([.03^2, .03^2, .03^2]);
Pv = diag([.01^2, .01^2, .01^2]);
Palpha = diag([0.01^2, 0.01^2, 0.01^2]);
Pba = diag([.01^2, 0.01^2, .01^2]); % just kind of making these up for now
Pbg = diag([.01^2, 0.01^2, .01^2]);
Ppf = diag([.01^2, .01^2, .01^2]);

%IMU
% see STIM300 ButterflyGyro Datasheet
g = 9.80665; % m/s^2
gyroParams = gyroparams('NoiseDensity', PSDg, 'BiasInstability', 0.3*pi/180/3600); % 0.3 deg/hr => rad/s
accelParams = accelparams('NoiseDensity', PSDa, 'BiasInstability', .05*.001*g); % taking unit mg to mean "milli-g's" => m/s^2
imu = imuSensor('SampleRate', 1/dt, 'Gyroscope', gyroParams, 'Accelerometer', accelParams);
[a_true_imu, w_true_imu] = I.trueInput(tVec);
% a_true_imu(3,:) = a_true_imu(3,:) - g; % looks like gravity is already inserted for you
[accelData, gyroData] = imu(a_true_imu', w_true_imu'); % probably have to move this into the loop once you start doing monte carlos
% plotIMU(accelData, gyroData, tVec, 2, 3)

%sensor
cameraRate = 5;
T_b2c = eye(3); %body to camera
Vc = diag([.001^2,.001^2]);

Nmonte = 1;
%%%Initialization
%measurement vectors
wmVec = zeros(T/dt,3,Nmonte);
amVec = zeros(T/dt,3,Nmonte);
xVec = zeros(T/dt,3,Nmonte);
vVec = zeros(T/dt,3,Nmonte);
qVec = zeros(T/dt,4,Nmonte);
baVec = zeros(T/dt,3,Nmonte);
bgVec = zeros(T/dt,3,Nmonte);

PxVec = zeros(T/dt,Nmonte);
PyVec = zeros(T/dt,Nmonte);
PzVec = zeros(T/dt,Nmonte);
PvxVec = zeros(T/dt,Nmonte);
PvyVec = zeros(T/dt,Nmonte);
PvzVec = zeros(T/dt,Nmonte);
PphiVec = zeros(T/dt,Nmonte);
PthVec = zeros(T/dt,Nmonte);
PpsiVec = zeros(T/dt,Nmonte);

x_err = zeros(T/dt,Nmonte);
y_err = zeros(T/dt,Nmonte);
z_err = zeros(T/dt,Nmonte);
vx_err = zeros(T/dt,Nmonte);
vy_err = zeros(T/dt,Nmonte);
vz_err = zeros(T/dt,Nmonte);
phiVec = zeros(T/dt,Nmonte);
thVec = zeros(T/dt,Nmonte);
psiVec = zeros(T/dt,Nmonte);

%%% OBSERVABILITY TESTING %%%
M = [];
nullSpaceDim = [];
Mrank = [];
PHI = eye(18);


for iMonte=1:Nmonte
    
%     a_true_i = 0.15*ones(3,1); % m/s^2
%     w_true_b = (pi/10)*ones(3,1); % rad/s
    [a0_true_i, w0_true_b] = I.trueInput(0);
%     T0_i2b = quat2dcm(q_init');
    
    wmVec(1,:,iMonte) = gyroData(1,:); %w0_true_b' + sqrt(PSDg/dt)*randn(3,1)';
    amVec(1,:,iMonte) = accelData(1,:); %a0_true_i' + sqrt(PSDa/dt)*randn(3,1)';

%     poseVec = []; % stores a [q,p] pair every time an image is recorded
    P = blkdiag(Px, Pv, Palpha, Pba, Pbg, Ppf);
    PxVec(1,iMonte) = P(1,1);
    PyVec(1,iMonte) = P(2,2);
    PzVec(1,iMonte) = P(3,3);
    PvxVec(1,iMonte) = P(4,4);
    PvyVec(1,iMonte) = P(5,5);
    PvzVec(1,iMonte) = P(6,6);
    PphiVec(1,iMonte) = P(7,7);
    PthVec(1,iMonte) = P(8,8);
    PpsiVec(1,iMonte) = P(9,9);
    
    
    x_i = x_init + chol(Px)*randn(3,1);  % inertial
    v_i = v_init + chol(Pv)*randn(3,1);  % inertial
    da = chol(Palpha)*randn(3,1);  % body
    nt = norm(da);
    if nt>1e-10
        dq = [cos(0.5*nt), ...
            (da/nt)'.*sin(0.5*nt)];
        q_i2b = quatmultiply(q_init',dq)';
    else
        q_i2b = q_init; % inertial-to-body
    end
    q_true_i2b = q_init;
    ba = chol(Pba)*randn(3,1); % biases assumed to be zero
    bg = chol(Pbg)*randn(3,1);
    pfVec_i = I.map + .01*randn(3,1); % add initial guess of single point location
    pfTags = [1];
    
    xVec(1,:,iMonte) = x_i';
    vVec(1,:,iMonte) = v_i';
    %aVec(1,:,iMonte) = a_init';
    qVec(1,:,iMonte) = q_i2b';
    x_err(1,iMonte) = x_i(1) - x_init(1);
    y_err(1,iMonte) = x_i(2) - x_init(2);
    z_err(1,iMonte) = x_i(3) - x_init(3);
    vx_err(1,iMonte) = v_i(1) - v_init(1);
    vy_err(1,iMonte) = v_i(2) - v_init(2);
    vz_err(1,iMonte) = v_i(3) - v_init(3);
    q_b2i = [q_i2b(1); -q_i2b(2:4)];
    dq_b = quatmultiply(q_b2i',q_init')';
    phiVec(1,iMonte) = 2*dq_b(2)/dq_b(1);
    thVec(1,iMonte) = 2*dq_b(3)/dq_b(1);
    psiVec(1,iMonte) = 2*dq_b(4)/dq_b(1);
%     V_Phi_test = zeros(3,T/dt);
%     X_Phi_test = zeros(3,T/dt);
    
    for ii = 2:(T/dt)
        t_curr = (ii-1)*dt;
        T_i2b_true = quat2dcm(q_true_i2b'); %inertial-to-body DCM
        
        %measurement
%         [a_true_i,w_true_b] = I.trueInput(t_curr);
%         a_true_b = T_i2b_true*a_true_i;
        am_b = T_i2b_true*accelData(ii,:)'; %T_i2b_true*a_true_i + sqrt(PSDa/dt)*randn(3,1); %body acceleration [m/s^2]
        wm_b = gyroData(ii,:)'; %w_true_b + sqrt(PSDg/dt)*randn(3,1); %body angular rate [rad/s]
        
        %propagation
        x_true = states_hist_true(ii,1:3)';
        vx_true = states_hist_true(ii,4:6)';
        q_true_i2b = states_hist_true(ii,7:10)';
        pf_true = I.map; % true point location
        [xP_i, vP_i, qP_i2b, pfVecP_i, PP, Phi] = prop(am_b, wm_b, x_true, vx_true, q_true_i2b, ba, bg, pf_true, dt, P, PSDa, PSDg, qf);
        %%%% PHI COMPARISON %%%%
%         g = [0;0;-1];
%         T_i2b_hat = quat2dcm(q_i2b'); %old inertial-to-body DCM
%         Phi_va = Phi(4:6,7:9);
%         Phi_xa = Phi(1:3,7:9);
%         V_Phi_test(:,ii) = (Phi_va*T_i2b_hat - cpe(v_i) + cpe(vP_i))*g;
%         X_Phi_test(:,ii) = (Phi_xa*T_i2b_hat - dt*cpe(v_i) - cpe(x_i) + cpe(xP_i))*g;

        %store propagated values and measurements
        amVec(ii,:,iMonte) = am_b';
        wmVec(ii,:,iMonte) = wm_b';
        xVec(ii,:,iMonte) = xP_i';
        vVec(ii,:,iMonte) = vP_i';
        qVec(ii,:,iMonte) = qP_i2b';
        
        %calculate position error and covariance in x, y, and z
        x_err(ii,iMonte) = xP_i(1) - x_true(1);
        y_err(ii,iMonte) = xP_i(2) - x_true(2);
        z_err(ii,iMonte) = xP_i(3) - x_true(3);
        PyVec(ii,iMonte) = PP(2,2);
        PzVec(ii,iMonte) = PP(3,3);
        
        %calculate velocity error and covariance: vx, vy, vz
        vx_err(ii,iMonte) = vP_i(1) - vx_true(1);
        vy_err(ii,iMonte) = vP_i(2) - vx_true(2);
        vz_err(ii,iMonte) = vP_i(3) - vx_true(3);
        PvxVec(ii,iMonte) = PP(4,4);
        PvyVec(ii,iMonte) = PP(5,5);
        PvzVec(ii,iMonte) = PP(6,6);
        
        %add an error and true quaternion
%         nt = norm(w_true_b)*dt;
%         if nt>1e-10
%             qdth_b = [cos(0.5*nt), ...
%                 (w_true_b*dt/nt)'.*sin(0.5*nt)];
%             q_true_i2b = quatmultiply(q_true_i2b',qdth_b)';
%         else
%             q_true_i2b = q_init;
%         end
        qP_b2i = [qP_i2b(1); -qP_i2b(2:4)];
        dq_b = quatmultiply(qP_b2i',q_true_i2b')';
        phiVec(ii,iMonte) = 2*dq_b(2)/dq_b(1);
        thVec(ii,iMonte) = 2*dq_b(3)/dq_b(1);
        psiVec(ii,iMonte) = 2*dq_b(4)/dq_b(1);
        PxVec(ii,iMonte) = PP(1,1);
        PphiVec(ii,iMonte) = PP(7,7);
        PthVec(ii,iMonte) = PP(8,8);
        PpsiVec(ii,iMonte) = PP(9,9);
        
        
        %update 
        if mod(ii,cameraRate) == 0
            [rm,pfTagsNew,pfTagsObserved,V] = rangeMeasurement(x_true, q_true_i2b, pfTags, I);
            nNewPts = length(pfTagsNew);
%             if nNewPts > 0
%                 Vf = .00000001;
%                 newPts = I.map(:,pfTagsNew) + sqrt(Vf)*randn(3,nNewPts); %noisy new features
%                 pfVec_i = [pfVec_i,newPts];
%                 pfTags = [pfTags,pfTagsNew];
%                 PP = blkdiag(PP,Vf*eye(3*nNewPts)); %if features were seen, augment covariance
%             end
            [xU_i, vU_i, qU_i2b, baU, bgU, pfVecU_i, PU, H] = update(rm(:), pfTagsObserved, pfTags, xP_i, vP_i, qP_i2b, ba, bg, ...
                pfVec_i, PP, V, I, x_true, q_true_i2b, pf_true); %%
            
            
            x_i = xU_i;
            v_i = vU_i;
            q_i2b = qU_i2b; 
            ba = baU; 
            bg = bgU;
            pfVec_i = pfVecU_i;
            P = PU;  
            
            % Rank test for observability matrix M
%             dimPHI = size(PHI,1);
%             dimH = size(H,2);
%             if (dimPHI - dimH) ~= 0
%                 PHI = blkdiag(PHI,eye(dimH-dimPHI));
%             end
            M_curr = H*PHI;
            M = [M; M_curr];
            Mrank = [Mrank; rank(M)];
            nullSpaceDim = [nullSpaceDim; size(M,2) - rank(M)];
            PHI = Phi*PHI;
        else
            x_i = xP_i;
            v_i = vP_i;
            q_i2b = qP_i2b; 
            pfVec_i = pfVecP_i;
            P = PP;  
        end
%         if ~isempty(pfVec_i)
%              figure(9); plot3(pfVec_i(1,:),pfVec_i(2,:),pfVec_i(3,:),'*'); axis equal; pause(.01)
%         end
    end
end

% Plotting
figure(4); clf;
subplot(3,1,1)
hold on
plot(dt:dt:T, x_err, 'b')
plot(dt:dt:T, sqrt(PxVec)*3, 'r')
plot(dt:dt:T, -sqrt(PxVec)*3, 'r')
hold off
title("x-position error and variance (3\sigma)")
ylabel("x [m]")
subplot(3,1,2)
hold on
plot(dt:dt:T, y_err, 'b')
plot(dt:dt:T, sqrt(PyVec)*3, 'r')
plot(dt:dt:T, -sqrt(PyVec)*3, 'r')
hold off
title("y-position error and variance (3\sigma)")
ylabel("y [m]")
subplot(3,1,3)
hold on
plot(dt:dt:T, z_err, 'b')
plot(dt:dt:T, sqrt(PzVec)*3, 'r')
plot(dt:dt:T, -sqrt(PzVec)*3, 'r')
hold off
title("z-position error and variance (3\sigma)")
ylabel("z [m]")

figure(5); clf;
subplot(3,1,1)
hold on
plot(dt:dt:T, vx_err, 'b')
plot(dt:dt:T, sqrt(PvxVec)*3, 'r')
plot(dt:dt:T, -sqrt(PvxVec)*3, 'r')
hold off
title("x-velocity error and variance (3\sigma)")
ylabel("v_x [m/s]")
subplot(3,1,2)
hold on
plot(dt:dt:T, vy_err, 'b')
plot(dt:dt:T, sqrt(PvyVec)*3, 'r')
plot(dt:dt:T, -sqrt(PvyVec)*3, 'r')
hold off
title("y-velocity error and variance (3\sigma)")
ylabel("v_y [m/s]")
subplot(3,1,3)
hold on
plot(dt:dt:T, vz_err, 'b')
plot(dt:dt:T, sqrt(PvzVec)*3, 'r')
plot(dt:dt:T, -sqrt(PvzVec)*3, 'r')
hold off
title("z-velocity error and variance (3\sigma)")
ylabel("v_z [m/s]")

figure(6); clf;
subplot(3,1,1)
hold on
plot(dt:dt:T, phiVec, 'b')
plot(dt:dt:T, sqrt(PphiVec)*3, 'r')
plot(dt:dt:T, -sqrt(PphiVec)*3, 'r')
hold off
title("\phi-attitude error and variance (3\sigma)")
ylabel("[rad]")
subplot(3,1,2)
hold on
plot(dt:dt:T, thVec, 'b')
plot(dt:dt:T, sqrt(PthVec)*3, 'r')
plot(dt:dt:T, -sqrt(PthVec)*3, 'r')
hold off
title("\theta-attitude error and variance (3\sigma)") 
ylabel("[rad]")
subplot(3,1,3)
hold on
plot(dt:dt:T, psiVec, 'b')
plot(dt:dt:T, sqrt(PpsiVec)*3, 'r')
plot(dt:dt:T, -sqrt(PpsiVec)*3, 'r')
hold off
title("\psi-attitude error and variance (3\sigma)")
ylabel("[rad]")

figure(1); %clf;
axisParam = I.plotAxisParam;
plot3(xVec(:,1), xVec(:,2), xVec(:,3), '-k')
hold on
for ii = 1:floor(length(tVec)/5):size(qVec,1)
    T_i2b = quat2dcm(qVec(ii,:,1));
    x = T_i2b'*[0.5*axisParam;0;0];
    y = T_i2b'*[0;0.5*axisParam;0];
    z = T_i2b'*[0;0;0.5*axisParam];
    quiver3(xVec(ii,1,1), xVec(ii,2,1), xVec(ii,3,1), x(1), x(2), x(3), 'r');
    quiver3(xVec(ii,1,1), xVec(ii,2,1), xVec(ii,3,1), y(1), y(2), y(3), 'g');
    quiver3(xVec(ii,1,1), xVec(ii,2,1), xVec(ii,3,1), z(1), z(2), z(3), 'b');
end
hold off
title("Trajectory of body frame (red = x_b, green = y_b, blue = z_b)")
axis equal

figure(10); plot(nullSpaceDim, 'LineWidth', 2); title('Dimension of nullspace of H-\Phi matrix')
figure(11); plot(Mrank, 'LineWidth', 2); title('Rank of M')

function [cm,pfTagsNew,pfTagsObserved,V] = cameraMeasurement(p,q,pfTags,I)
%CAMERAMEASUREMENT Generates a simulated camera measurement of features in
%I.map
Vc = .001;
T_i2b = quat2dcm(q'); 
pf_mat = I.map; %feature positions
nPts = I.nPts;
p_mat = repmat(p,1,nPts);
map_b = T_i2b*(pf_mat - p_mat);
mask = (map_b(3,:) > 0) & (abs(map_b(2,:)./map_b(3,:)) <= 1.75) ...
    & (abs(map_b(1,:)./map_b(3,:)) <= 1.75); %points in front of the camera, roughly 60 deg FOV up/down left/right
map_b = map_b(:,mask); 
nPtsInView = size(map_b,2);
cm = map_b(1:2,:)./map_b(3,:) + Vc*randn(2,nPtsInView); %add noise
pfTagsObserved = find(mask);
pfTagsNew = setdiff(pfTagsObserved,pfTags); % Did we observe any new features?
V = Vc^2*eye(2*nPtsInView);
% plotCam(cm,mask,7);
end

function [zHat,H] = cameraMeasModel(p,q,pfVec,pfTags,pfTagsObserved,I)
%CAMERAMEASMODEL Generates H-matrix for camera measurement based on
%current estimated state of the robot
T_i2b = quat2dcm(q'); 
% pf_mat = I.map; %feature positions
nPtsInState = length(pfTags);
nPtsInView = length(pfTagsObserved);
indObserved = zeros(1,nPtsInView);
for ii = 1:nPtsInView
    indObserved(ii) = find(pfTags==pfTagsObserved(ii));
end
pfVecObserved = pfVec(:,indObserved);
p_mat = repmat(p,1,nPtsInView);
pfVec_b = T_i2b*(pfVecObserved - p_mat);
% pfVec_b = pfVec_b(:,mask);
% nPtsInView = size(pfVec_b,2);
H = zeros(2*nPtsInView,15+3*nPtsInState);
for ii = 1:nPtsInView
    px = pfVec_b(1,ii); py = pfVec_b(2,ii); pz = pfVec_b(3,ii);
    Hc = (1/pz^2)*[pz, 0, -px; 0, pz, -py]; 
    Hth = cpe([px;py;pz]);
    Hp = -T_i2b;
    Hf = T_i2b;
    H(2*ii-1:2*ii,1:9) = Hc*[Hp, zeros(3), Hth];
    ind = 3*indObserved(ii)-2:3*indObserved(ii);
    H(2*ii-1:2*ii,15+ind) = Hc*Hf;
end
zHat = pfVec_b(1:2,:)./pfVec_b(3,:);   
zHat = zHat(:);
end

function [rm,pfTagsNew,pfTagsObserved,V] = rangeMeasurement(p,q,pfTags,I)
%RANGEMEASUREMENT Generates a simulated range measurement of features in
%I.map
Vr = .01; % 1 cm of range noise
T_i2b = quat2dcm(q'); 
pf_mat = I.map; %feature positions
nPts = I.nPts;
p_mat = repmat(p,1,nPts);
map_b = T_i2b*(pf_mat - p_mat);
mask = [1]; % only looking at 1 feature
%     (map_b(3,:) > 0) & (abs(map_b(2,:)./map_b(3,:)) <= 1.75) ...
%     & (abs(map_b(1,:)./map_b(3,:)) <= 1.75); %points in front of the range measurement device, roughly 60 deg FOV up/down left/right
map_b = map_b(:,mask); 
nPtsInView = size(map_b,2);
rm = sqrt(map_b(1,:).^2 + map_b(2,:).^2 + map_b(3,:).^2) + Vr*randn();
pfTagsObserved = find(mask);
pfTagsNew = setdiff(pfTagsObserved,pfTags); % Did we observe any new features?
V = Vr^2*eye(nPtsInView);
end

function [zHat,H] = rangeMeasModel(p,q,pfVec,pfTags,pfTagsObserved,I)
%CAMERAMEASMODEL Generates H-matrix for camera measurement 
% T_i2b = quat2dcm(q'); 
% pf_mat = I.map; %feature positions
nPtsInState = length(pfTags);
nPtsInView = length(pfTagsObserved);
indObserved = zeros(1,nPtsInView);
for ii = 1:nPtsInView
    indObserved(ii) = find(pfTags==pfTagsObserved(ii));
end
pfVecObserved = pfVec(:,indObserved);
% p_mat = repmat(p,1,nPtsInView);
% pfVec_b = T_i2b*(pfVecObserved - p_mat);
% pfVec_b = pfVec_b(:,mask);
% nPtsInView = size(pfVec_b,2);
H = zeros(nPtsInView,15+3*nPtsInState);
for ii = 1:nPtsInView
    pf = pfVecObserved(:,ii);
    hr = 1/sqrt((pf-p)'*(pf-p));
    H(ii,1:3) = -hr*(pf-p)';
    ind = 3*indObserved(ii)-2:3*indObserved(ii);
    H(ii,15+ind) = hr*(pf-p)';
end
pfVecDiff = pfVecObserved - p;
zHat = sqrt(pfVecDiff(1,:).^2 + pfVecDiff(2,:).^2 + pfVecDiff(3,:).^2);
zHat = zHat(:);
end

function [xU_i, vU_i, qU_i2b, baU, bgU, pfVecU_i, PU, H] = update(z, pfTagsObserved, pfTags, x_i, v_i, q_i2b, ba, bg, pfVec_i,P,V,I,x_true,q_true_i2b,pf_true)
[zHat,H] = rangeMeasModel(x_true,q_true_i2b,pf_true,pfTags,pfTagsObserved,I);
K = (P*H')/(H*P*H' + V);
delStates = K*(z-zHat);
xU_i = x_i + delStates(1:3);
vU_i = v_i + delStates(4:6);
delStates(10:15)'
baU = ba + delStates(10:12);
bgU = bg + delStates(13:15);
delpf = reshape(delStates(16:end),3,[]);
pfVecU_i = pfVec_i;
nPtsInView = length(pfTagsObserved); %only update the features that were observed
for ii = 1:nPtsInView
    indObserved = find(pfTags==pfTagsObserved(ii));
    pfVecU_i(:,indObserved) = pfVecU_i(:,indObserved) + delpf(:,ii);
end
%%%%% ADDITIVE UPDATE %%%%%
% qU_i2b = q_i2b + [0;delStates(7:9)];
%%%%%%
%%%%%% MULITPLICATIVE UPDATE %%%%%
dth = delStates(7:9);
nt = norm(dth);
qdth = [cos(0.5*nt), (dth/nt)'.*sin(0.5*nt)];
qU_i2b = quatmultiply(q_i2b',qdth)'; %quaternion update (Sola 64)
%%%%%
PU = P - K*(H*P*H'+V)*K'; %(eye(9) - K*H)*P*(eye(9) - K*H)' + K*V*K'; %Joseph form covariance update (Sola 63)
end


function [xP_i, vP_i, qP_i2b, pfVecP_i, PP, Phi] = prop(am_b, wm_b, x_i, v_i, q_i2b, ba, bg, pfVec_i, dt, P, PSDa, PSDw, qf)
% x = [x_i, v_i, a_ib]
% x_i : 3x1 IMU position in the inertial frame
% v_i : 3x1 IMU velocity in the inertial frame
% a_b : 3x1 estimated error angle (since last measurement update)
% TBIhat : 3x3 estimated DCM from body-frame to inertial frame
% am_b : measured acceleration in body frame
% wm_b : measured angular rate of the body frame
% G = 3x3 gravity matrix
%
% Covariance propagation:
% P(k+1) = Phi*P*Phi' + B*Q*B'

% Pre-processing
% G = zeros(3,3);
g = 9.80665; % m/s^2
T_i2b_hat = quat2dcm(q_i2b'); %inertial-to-body DCM
am_i = T_i2b_hat' * (am_b - ba);
dth_b = (wm_b - bg) * dt; %delta theta [rad]

% Propagate states
xP_i = x_i + v_i*dt + 0.5*am_i*dt^2;
vP_i = v_i + am_i*dt;
nt = norm(dth_b);
if nt>1e-10
    qdth_b = [cos(0.5*nt), ...
        (dth_b/nt)'.*sin(0.5*nt)];
    qP_i2b = quatmultiply(q_i2b',qdth_b)';
else
    qP_i2b = q_i2b;
end
pfVecP_i = pfVec_i;
nPts = size(pfVec_i,2);

% Uncertainties
Qg = dt*PSDw*eye(3);
Qbg = dt*(0.3*pi/180/3600)*eye(3); % copying these numbers in from the imu for now
Qa = dt*PSDa*eye(3); 
Qba = dt*(.05*.001*g)*eye(3);% copying these numbers in from the imu for now
Q = blkdiag(Qg, Qbg, Qa, Qba);

% Propagate covariance
% See Hesch et al 2014
F_x_x = zeros(3); F_x_v = eye(3); F_x_a = zeros(3); F_x_ba = zeros(3); F_x_bg = zeros(3); 
F_v_x = zeros(3); F_v_v = zeros(3); F_v_a = -T_i2b_hat'*cpe(am_i); F_v_ba = -T_i2b_hat'; F_v_bg = zeros(3);
F_a_x = zeros(3); F_a_v = zeros(3); F_a_a = -cpe(wm_b - bg); F_a_ba = zeros(3); F_a_bg = -eye(3);
F_ba_x = zeros(3); F_ba_v = zeros(3); F_ba_a = zeros(3); F_ba_ba = zeros(3); F_ba_bg = zeros(3); 
F_bg_x = zeros(3); F_bg_v = zeros(3); F_bg_a = zeros(3); F_bg_ba = zeros(3); F_bg_bg = zeros(3); 
F = [F_x_x, F_x_v, F_x_a, F_x_ba, F_x_bg; 
     F_v_x, F_v_v, F_v_a, F_v_ba, F_v_bg; 
     F_a_x, F_a_v, F_a_a, F_a_ba, F_a_bg;
     F_ba_x, F_ba_v, F_ba_a, F_ba_ba, F_ba_bg; 
     F_bg_x, F_bg_v, F_bg_a, F_bg_ba, F_bg_bg];
Phi = expm(F*dt)*eye(15); %xdot = Ax => x = e^(At)x0.
if nPts > 0
    Phi = blkdiag(Phi,eye(3*nPts));
end

G_Qg_x = zeros(3); G_Qbg_x = zeros(3); G_Qa_x = zeros(3); G_Qba_x = zeros(3); 
G_Qg_v = zeros(3); G_Qbg_v = zeros(3); G_Qa_v = -T_i2b_hat'; G_Qba_v = zeros(3);
G_Qg_a = -eye(3); G_Qbg_a = zeros(3);  G_Qa_a = zeros(3); G_Qba_a = zeros(3); 
G_Qg_ba = zeros(3); G_Qbg_ba = zeros(3); G_Qa_ba = zeros(3);  G_Qba_ba = eye(3);
G_Qg_bg = zeros(3); G_Qbg_bg = eye(3); G_Qa_bg = zeros(3); G_Qba_bg = zeros(3); 
 
G = [G_Qg_x, G_Qbg_x, G_Qa_x, G_Qba_x;
    G_Qg_v, G_Qbg_v, G_Qa_v, G_Qba_v;
    G_Qg_a, G_Qbg_a,  G_Qa_a, G_Qba_a; 
    G_Qg_ba, G_Qbg_ba, G_Qa_ba, G_Qba_ba;
    G_Qg_bg, G_Qbg_bg, G_Qa_bg, G_Qba_bg;
    zeros(3*nPts,12)];
     
% B_xna = -0.5*dt*T_i2b_hat';
% B_xng = zeros(3,3);
% B_vna = -T_i2b_hat';
% B_vng = zeros(3,3);
% B_ana = zeros(3,3);
% B_ang = eye(3,3);
% B = [B_xna, B_xng;...
%     B_vna, B_vng; ...
%     B_ana, B_ang; ...
%     zeros(3*nPts,6)];

PP = Phi*P*Phi' + G*Q*G';

end

function [statesDot] = trueKinematics(t,states,I,qf)
[a_i,w_b] = I.trueInput(t);
v = states(4:6);
q = states(7:10);
xdot = v;
vdot = a_i;
% qdot = 0.5*qf.Omega(w_b)*q; 
qdot = quatmultiply(q',[0,0.5*w_b'])';
statesDot = [xdot;vdot;qdot];
end

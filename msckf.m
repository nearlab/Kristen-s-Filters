%%%Inputs
close all; clear
I = generateHeschInput();
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
gVec_i = [0;0;-9.81];

states_init = [x_init; v_init; q_init];
[t_hist,states_hist_true] = ode45(@(t,states)trueKinematics(t,states,I,qf),tVec,states_init);
plotTrajectory(t_hist,states_hist_true,1,I)
%uncertainties
%x: +/- .01 m
%v: +/- .01 m/s
%q: +/- .1 rad
PSDa = (.7/60)^2; % [(m/s)^2/s] = 0.7 m/s/sqrt(hr) %%
PSDg = (0.15*pi/180/60)^2; % [rad^2/s] = 0.15 deg/sqrt(hr) %%
Pmap = diag([.03^2, .03^2, .03^2]);
Px = diag([.03^2, .03^2, .03^2]);
Pv = diag([.01^2, .01^2, .01^2]);
Palpha = diag([0.01^2, 0.01^2, 0.01^2]);

%sensor
cameraRate = 10; %%
T_b2c = eye(3); %body to camera 
q_b2c = [1,0,0,0]'; 
x_cb = [0,0,0]'; %position of camera in body frame
sig_Im = .01;
Vc = diag([sig_Im^2,sig_Im^2]);

Nmonte = 5;
%%%Initialization
%measurement vectors
wmVec = zeros(T/dt,3,Nmonte);
amVec = zeros(T/dt,3,Nmonte);
xVec = zeros(T/dt,3,Nmonte);
vVec = zeros(T/dt,3,Nmonte);
qVec = zeros(T/dt,4,Nmonte);

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

for iMonte=1:Nmonte
    M = msckfDataManager(); % re-initialize the data manager on each run
    map_error = chol(Pmap)*randn(3,1);
    q_true_i2b = states_hist_true(1,7:10)';
    T_i2b_true = quat2dcm(q_true_i2b');
    [a0_true_i, w0_true_b] = I.trueInput(0);
    
    wmVec(1,:,iMonte) = w0_true_b' + sqrt(PSDg/dt)*randn(3,1)'; 
    amVec(1,:,iMonte) = (T_i2b_true*(a0_true_i - gVec_i))' + sqrt(PSDa/dt)*randn(3,1)'; 

    P = blkdiag(Px + Pmap, Pv, Palpha);
    PxVec(1,iMonte) = P(1,1);
    PyVec(1,iMonte) = P(2,2);
    PzVec(1,iMonte) = P(3,3);
    PvxVec(1,iMonte) = P(4,4);
    PvyVec(1,iMonte) = P(5,5);
    PvzVec(1,iMonte) = P(6,6);
    PphiVec(1,iMonte) = P(7,7);
    PthVec(1,iMonte) = P(8,8);
    PpsiVec(1,iMonte) = P(9,9);
    
    
    x_i = x_init + chol(Px)*randn(3,1) + map_error;  % inertial %%
    v_i = v_init + chol(Pv)*randn(3,1);  % inertial %%
    da = chol(Palpha)*randn(3,1);  % body %%
    nt = norm(da);
    if nt>1e-10
        dq = [cos(0.5*nt), ...
            (da/nt)'.*sin(0.5*nt)];
        q_i2b = quatnormalize(quatmultiply(q_init',dq))';
    else
        q_i2b = q_init; % inertial-to-body
    end
    q_true_i2b = q_init;
    poseVec_i2c = []; % stores a [q,p] pair every time an image is recorded
%     pfTags = [];
    updateFlag = 0;
    
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
        t_curr = (ii-1)*dt; % propagating from time (ii-1)*dt to time (ii)*dt
%         x_true_i = states_hist_true(ii,1:3)';
        T_i2b_true = quat2dcm(states_hist_true(ii-1,7:10)); %inertial-to-body DCM
        
        %measurement
        [a_true_i,w_true_b] = I.trueInput((ii-2)*dt);
%         a_true_b = T_i2b_true*a_true_i;
        am_b = T_i2b_true*(a_true_i - gVec_i) + sqrt(PSDa/dt)*randn(3,1); %body acceleration [m/s^2] 3
        wm_b = w_true_b + sqrt(PSDg/dt)*randn(3,1); %body angular rate [rad/s] 
        
        %propagation
        [xP_i, vP_i, qP_i2b, poseVecP_i2c, PP] = prop(am_b, wm_b, x_i, v_i, q_i2b, poseVec_i2c, dt, P, PSDa, PSDg, qf);
%         am_b = T_i2b_true*(vP_i-v_i)/dt;
        
%         L = eig(PP);
%         if ~isreal(L) || min(real(L)) < 0
%             L;
%         end
%         T_i2b_true = quat2dcm(q_true_i2b');
        
        %store propagated values and measurements
        amVec(ii,:,iMonte) = am_b';
        wmVec(ii,:,iMonte) = wm_b';
        xVec(ii,:,iMonte) = xP_i';
        vVec(ii,:,iMonte) = vP_i';
        qVec(ii,:,iMonte) = qP_i2b';
        
        %calculate position error and covariance in x, y, and z
        x_true = states_hist_true(ii,1:3)';
        x_err(ii,iMonte) = xP_i(1) - x_true(1);
        y_err(ii,iMonte) = xP_i(2) - x_true(2);
        z_err(ii,iMonte) = xP_i(3) - x_true(3);
        PxVec(ii,iMonte) = PP(1,1);
        PyVec(ii,iMonte) = PP(2,2);
        PzVec(ii,iMonte) = PP(3,3);
        
        %calculate velocity error and covariance: vx, vy, vz
        vx_true = states_hist_true(ii,4:6)';
        vx_err(ii,iMonte) = vP_i(1) - vx_true(1);
        vy_err(ii,iMonte) = vP_i(2) - vx_true(2);
        vz_err(ii,iMonte) = vP_i(3) - vx_true(3);
        PvxVec(ii,iMonte) = PP(4,4);
        PvyVec(ii,iMonte) = PP(5,5);
        PvzVec(ii,iMonte) = PP(6,6);
        
        %add an error and true quaternion
%         q_true_i2b = states_hist_true(ii,7:10)';
%         nt = norm(w_true_b)*dt;
%         if nt>1e-10
%             qdth_b = [cos(0.5*nt), ...
%                 (w_true_b*dt/nt)'.*sin(0.5*nt)];
%             q_true_i2b = quatmultiply(q_true_i2b',qdth_b)';
%         else
%             q_true_i2b = q_init;
%         end
        q_true_i2b = states_hist_true(ii,7:10)';
        qP_b2i = [qP_i2b(1); -qP_i2b(2:4)];
        dq_b = quatmultiply(qP_b2i',q_true_i2b')';
        phiVec(ii,iMonte) = 2*dq_b(2)/dq_b(1);
        thVec(ii,iMonte) = 2*dq_b(3)/dq_b(1);
        psiVec(ii,iMonte) = 2*dq_b(4)/dq_b(1);
        PphiVec(ii,iMonte) = PP(7,7);
        PthVec(ii,iMonte) = PP(8,8);
        PpsiVec(ii,iMonte) = PP(9,9);
        
        
        %update 
        if mod(ii,cameraRate) == 0
%             % Store pose
%             T_i2b = quat2dcm(qP_i2b');
%             q_i2c = quatmultiply(qP_i2b',q_b2c')';
%             x_ci = xP_i + T_i2b'*x_cb; %position of camera in inertial frame
%             pose = [q_i2c; x_ci]; %length 7 addition to state vector
%             
%             %%%%% RZ messed this up
% %             x_ci = x_true_i + T_i2b_true'*x_cb; %position of camera in inertial frame
% %             pose = [q_true_i2b; x_ci]; %length 7 addition to state vector
%             %%%%%%
%             
%             
%             poseVecP_i2c = [poseVecP_i2c, pose]; % append new pose to the state vector each time an image is taken
%             M.poseTimestampList = [M.poseTimestampList, t_curr]; % each pose associated with unique timestamp
%             M.nPoses = M.nPoses + 1;
%             
%             % Augment covariance matrix to include new pose
%             N = M.nPoses;
%             T_i2c = quat2dcm(q_i2c');
%             J = [zeros(3,6), T_b2c, zeros(3,6*(N-1));
%                  eye(3), zeros(3), qf.cpe(T_i2c'*x_ci), zeros(3,6*(N-1))];
%             PP = [eye(6*(N-1)+9); J]*PP*[eye(6*(N-1)+9); J]';
% %             PP = sqrtm(PP*PP);
% %             [R,flag] = chol(PP);
%             L = eig(PP);
%             if ~isreal(L) || min(real(L)) < 0
%                 t_curr
%                 min(real(L))
% %                 eig(PP)
%                 L;
%             end
%             % If there are more than 7 poses, throw out the oldest one
%             if M.nPoses > 7
%                  poseVecP_i2c = poseVecP_i2c(:,2:end); % throw out the oldest pose
%                  M.nPoses = M.nPoses - 1; 
%                  M.poseTimestampList(1) = [];
%                  PP = PP([1:9,16:end],[1:9,16:end]);
%             end
            
            % Measure
            [cm,pfTagsNew,pfTagsObserved,V] = cameraMeasurement(x_true, q_true_i2b, M.tagList, I);
            nNewPts = length(pfTagsNew);
            if nNewPts > 0
                newPtPkgs = []; %weird to initialize struct arrays
                for jj = 1:nNewPts
                    cmInd = find(pfTagsObserved == pfTagsNew(jj));
                    newPtPkg = pointPackage(I,pfTagsNew(jj),cm(:,cmInd),t_curr); %, pose);
                    newPtPkgs = [newPtPkgs, newPtPkg];
                end
                M.pfPkgList = [M.pfPkgList, newPtPkgs];
                M.tagList = [M.tagList, pfTagsNew];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Initialize the state perfectly and set noises to zero
            % Use your estimated pose and xIMU to compute the estimated
            % measurements
            % Make sure the actual measurements match the estimated ones
            % (ie turn off all the updates, all the noises; make sure
            % estimated measurements match real measurements)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Add measurement to points seen before
            pfTagsOld = setdiff(pfTagsObserved,pfTagsNew);
            nOldPts = length(pfTagsOld);
            if nOldPts > 0
                for kk = 1:nOldPts
                    pkgInd = find(M.tagList == pfTagsOld(kk));
                    cmInd = find(pfTagsObserved == pfTagsOld(kk));
                    currPkg = M.pfPkgList(pkgInd);
                    currPkg.t_list = [currPkg.t_list, t_curr];
                    currPkg.meas_list = [currPkg.meas_list, cm(:,cmInd)]; 
%                     currPkg.pose_list = [currPkg.pose_list, pose];
                    M.pfPkgList(pkgInd) = currPkg;
                end
            end
            
            % When a feature that has been tracked in a number of images is
            % no longer detected, then all the measurements of that feature
            % are processed at the same time (see MSCKF 2007).
            pfTagsOOF = setdiff(M.tagList,pfTagsObserved); %OOF = "out of frame"
            nPfTagsOOF = length(pfTagsOOF);
            if nPfTagsOOF > 0
                updateFlag = 1;
                for mm = 1:nPfTagsOOF 
                    pkgInd = find(M.tagList == pfTagsOOF(mm));
                    currPkg = M.pfPkgList(pkgInd);
                    % Least-squares estimate feature location based on
                    % observation data in the point package
%                     lsEstimatePf(currPkg,poseVec_i2c,M)
                    pfHat_i = currPkg.p_true + sqrt(1e-6)*randn(3,1) + map_error; 
%                     L = eig(PP);
%                     if ~isreal(L) || min(real(L)) < 0
%                         L;
%                     end
                    if mm == 1
                        [xU_i, vU_i, qU_i2b, poseVecU_i2c, PU] = update(currPkg, pfHat_i, xP_i, vP_i, qP_i2b, poseVecP_i2c, PP, sig_Im, I, M, qf); 
                    else
                        [xU_i, vU_i, qU_i2b, poseVecU_i2c, PU] = update(currPkg, pfHat_i, xU_i, vU_i, qU_i2b, poseVecU_i2c, PU, sig_Im, I, M, qf); 
                    end
                    % Remove the feature package from M
                    M.tagList(pkgInd) = [];
                    M.pfPkgList(pkgInd) = [];
                end
            end
            
            if updateFlag == 1
                x_i = xU_i;
                v_i = vU_i;
                q_i2b = qU_i2b; 
                poseVec_i2c = poseVecU_i2c;
                P = PU;  
            else 
                x_i = xP_i;
                v_i = vP_i;
                q_i2b = qP_i2b; 
                poseVec_i2c = poseVecP_i2c;
                P = PP;  
            end
            
            updateFlag = 0;
            
            % Store pose
            T_i2b = quat2dcm(q_i2b');
            q_i2c = quatmultiply(q_i2b',q_b2c')';
            x_ci = x_i + T_i2b'*x_cb; %position of camera in inertial frame
            pose = [q_i2c; x_ci]; %length 7 addition to state vector
            
            %%%%% RZ messed this up
%             x_ci = x_true_i + T_i2b_true'*x_cb; %position of camera in inertial frame
%             pose = [q_true_i2b; x_ci]; %length 7 addition to state vector
            %%%%%%
            
            poseVec_i2c = [poseVec_i2c, pose]; % append new pose to the state vector each time an image is taken
            M.poseTimestampList = [M.poseTimestampList, t_curr]; % each pose associated with unique timestamp
            M.nPoses = M.nPoses + 1;
            
            % Augment covariance matrix to include new pose
            N = M.nPoses;
            T_i2c = quat2dcm(q_i2c');
            J = [zeros(3,6), T_b2c, zeros(3,6*(N-1));
                 eye(3), zeros(3), cpe(T_i2b'*x_cb), zeros(3,6*(N-1))];  %%%
            P = [eye(6*(N-1)+9); J]*P*[eye(6*(N-1)+9); J]';
%             P = sqrtm(P*P');
%             [R,flag] = chol(PP);
            L = eig(P);
            if ~isreal(L) || min(real(L)) < 0
%                 t_curr
%                 min(real(L))
%                 eig(PP)
                L;
            end
            % If there are more than 7 poses, throw out the oldest one
            if M.nPoses > 20
                 poseVec_i2c = poseVec_i2c(:,2:end); % throw out the oldest pose
                 M.nPoses = M.nPoses - 1; 
                 M.poseTimestampList(1) = [];
                 P = P([1:9,16:end],[1:9,16:end]);
            end
            
        else
            x_i = xP_i;
            v_i = vP_i;
            q_i2b = qP_i2b; 
            poseVec_i2c = poseVecP_i2c;
            P = PP;  
        end
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
for ii = 1:floor(length(tVec)/20):size(qVec,1)
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


function [cm,pfTagsNew,pfTagsObserved,V] = cameraMeasurement(p,q,pfTags,I)
%CAMERAMEASUREMENT Generates a simulated camera measurement of features in
%I.map
Vc = .01;
T_i2b = quat2dcm(q'); 
pf_mat = I.map; %feature positions
nPts = I.nPts;
% p_mat = repmat(p,1,nPts);
map_b = T_i2b*(pf_mat - p);
fovParam = 0.5;
mask = (map_b(3,:) > 0) & (abs(map_b(2,:)./map_b(3,:)) <= fovParam) ...
    & (abs(map_b(1,:)./map_b(3,:)) <= fovParam); %points in front of the camera, roughly 45 deg FOV up/down left/right
map_b = map_b(:,mask); 
nPtsInView = size(map_b,2);
cm = map_b(1:2,:)./map_b(3,:)+ Vc*randn(2,nPtsInView); %add noise 
% cm(1,:) = -cm(1,:); %%%
pfTagsObserved = find(mask);
pfTagsNew = setdiff(pfTagsObserved,pfTags); % Did we observe any new features?
V = Vc^2*eye(2*nPtsInView);
% plotCam(cm,pfTagsObserved,fovParam,7);
end

function [r0,H0] = cameraMeasModel(ptPkg,pfHat_i,poseVec_i2c,qf,M)
%CAMERAMEASMODEL Generates H-matrices for point packages based on estimated
%states at times feature was observed
t_list = ptPkg.t_list;
nt = length(t_list);
pose_inds = [];
meas_inds = [];
for ii = 1:nt
    t = t_list(ii);
    ind = find(M.poseTimestampList == t);
    if ~isempty(ind)
        pose_inds = [pose_inds, ind];
        meas_inds = [meas_inds, ii];
    end
end
N = size(poseVec_i2c,2);
n = length(pose_inds);

% zHat = zeros(2,n);zHat = zeros(2,n); % inefficient! remove later
% for ii = 1:n
%     q = poseVec_i2c(1:4,pose_inds(ii));
%     p = poseVec_i2c(5:7,pose_inds(ii));
%     T_i2c = quat2dcm(q'); 
%     pfHat_c = T_i2c*(pfHat_i - p);
%     zHat(:,ii) = (1/pfHat_c(3))*[pfHat_c(1); pfHat_c(2)];
% end
% zHat = zHat(:);
% plotMeasModel(zHat,ptPkg.tag,7)

% See MSCKF Tech. Rpt. 2006 for Hx, Hf equations
Hx = zeros(2*n, 9+6*N);
Hf = zeros(2*n, 3);
zHat = zeros(2,n);
for ii = 1:n
    q = poseVec_i2c(1:4,pose_inds(ii));
    p = poseVec_i2c(5:7,pose_inds(ii)); %%%
    T_i2c = quat2dcm(q');
    pfHat_c = T_i2c*(pfHat_i - p);
    zHat(:,ii) = (1/pfHat_c(3))*[pfHat_c(1); pfHat_c(2)];
    J = 1/pfHat_c(3)*[1, 0, -pfHat_c(1)/pfHat_c(3); ...
                      0, 1, -pfHat_c(2)/pfHat_c(3)];
    idx = 9+(6*pose_inds(ii)-5:6*pose_inds(ii));
    Hx(2*ii-1:2*ii,idx) = [J*qf.cpe(pfHat_c), -J*T_i2c];
    Hf(2*ii-1:2*ii,:) = J*T_i2c;
end
zHat = zHat(:);
% plotMeasModel(zHat,ptPkg.tag,7)

% Find relevant measurements
z = ptPkg.meas_list(:,meas_inds);
z = z(:);

% Form matrix H0 and residual r0
r = z-zHat;
% figure(100); hold on; plot(t_list(end), r, '*');
% Givens rotations formulation
H0 = zeros(size(Hx));
r0 = zeros(size(r));
for ii = 1:3
    for jj = size(Hf,1):-1:2
        [G,~] = planerot(Hf(jj-1:jj,ii));
        Hf(jj-1:jj,ii) = G*Hf(jj-1:jj,ii);
        H0(jj-1:jj,:) = G*Hx(jj-1:jj,:);
        r0(jj-1:jj) = G*r(jj-1:jj);
    end
end
% H0 = H0(4:end-3,:);
% r0 = r0(4:end-3,:);
% Nullspace formulation
% A = null(Hf'); % left nullspace (i.e. Hf'*A = 0 ==> A'*Hf = 0)
% H0 = A'*Hx;
% r0 = A'*Hf;
% No adjustment for state correlation to errors in pf
% H0 = Hx;
% r0 = r;
end

function [xU_i, vU_i, qU_i2b, poseVecU_i2c, PU] = update(currPkg, pfHat_i, x_i, v_i, q_i2b, poseVec_i2c,P,sig_Im,I,M,qf)
% Find poses in poseVec associated with currPkg
t_list = currPkg.t_list; 
n = length(t_list); % number of poses associated with feature in ptPkg (= number of
% measurements in the ptPkg)
N = size(poseVec_i2c,2);

[r0,H0] = cameraMeasModel(currPkg, pfHat_i, poseVec_i2c, qf, M);
% z = currPkg.meas_list(:);
V = (sig_Im^2)*eye(size(r0,1)); % (sig_Im^2)*eye(size(A,2)); % 
K = (P*H0')/(H0*P*H0' + V);
% K(10:end,:)=0;
delStates = K*r0;
xU_i = x_i + delStates(1:3);
vU_i = v_i + delStates(4:6);
dth = delStates(7:9);
nt = norm(dth);
if nt>1e-10
    qdth = [cos(0.5*nt), (dth/nt)'.*sin(0.5*nt)];
else
    qdth = [1 0 0 0];
end
if nt > 0.25 || any(abs(delStates(1:6)) > 1)
    delStates;
end
qU_i2b = quatnormalize(quatmultiply(q_i2b',qdth))'; %quaternion update (Sola 64)

delPoseVec = reshape(delStates(10:end),6,[]);
poseVecU_i2c = poseVec_i2c;
for ii = 1:M.nPoses
    dth = delPoseVec(1:3,ii);
    nt = norm(dth);
    if nt>1e-10
        qdth = [cos(0.5*nt), (dth/nt)'.*sin(0.5*nt)];
    else
        qdth = [1 0 0 0];
    end
    poseVecU_i2c(1:4,ii) = quatnormalize(quatmultiply(poseVec_i2c(1:4,ii)',qdth))';
%     if norm(poseVecU_i2c(1:4,ii)) ~= 1
%         norm(poseVecU_i2c(1:4,ii))
%     end
    poseVecU_i2c(5:7,ii) = poseVec_i2c(5:7,ii) + delPoseVec(4:6,ii);
end
N = 9 + 6*size(poseVec_i2c,2);
PU = (eye(N) - K*H0)*P*(eye(N) - K*H0)' + K*V*K'; %Joseph form covariance update (Sola 63)
% L = eig(PU);
% if ~isreal(L) || min(real(L)) < 0
%     L;
% end
% Non-update!
% xU_i = x_i;
% vU_i = v_i;
% qU_i2b = q_i2b;
% poseVecU_i2c = poseVec_i2c;
% PU = P;
end


function [xP_i, vP_i, qP_i2b, poseVecP_i2c, PP] = prop(am_b, wm_b, x_i, v_i, q_i2b, poseVec_i2c, dt, P, PSDa, PSDw, qf)
% x = [x_i, v_i, a_ib]
% x_i : 3x1 IMU position in the inertial frame
% v_i : 3x1 IMU velocity in the inertial frame
% a_b : 3x1 estimated error angle (since last measurement update)
% TBIhat : 3x3 estimated DCM from body-frame to inertial frame
% am_b : measured acceleration in body frame
% wm_b : measured angular rate of the body frame
% G = 3x3 gravity matrix
%
% Phi = [Phi_xx Phi_vx Phi_ax
%        Phi_vx Phi_vv Phi_av
%        Phi_ax Phi_av Phi_aa]
% Phi_xx : I3x3 + (1/2)*(dt^2)*G
% Phi_vx : dt*I
% Phi_ax : -(1/2)*(dt^2)*TBIhat*[am_b x]
% Phi_vx : dt*G
% Phi_vv : I3x3
% Phi_av : -dt*TBIhat*[am_b x]
% Phi_ax : 03x3
% Phi_av : 03x3
% Phi_aa : I - cpe(wm_b)
%
% B = [B_xx B_xv
%      B_vx Bvv];
% B_xx = -(1/2)*(dt^2)*TBIhat;
% B_xv = 03x3
% B_vx = 03x3
% B_vv = I3x3
%
% Covariance propagation:
% P(k+1) = Phi*P*Phi' + B*Q*B'
gVec_i = [0;0;-9.81];

% Pre-processing
T_i2b_hat = quat2dcm(q_i2b'); %inertial-to-body DCM
am_i = T_i2b_hat' * am_b + gVec_i;
dth_b = wm_b * dt; %delta theta [rad]

% Propagate states
xP_i = x_i + v_i*dt + 0.5*am_i*dt^2;
vP_i = v_i + am_i*dt;
nt = norm(dth_b);
if nt>1e-10
    qdth_b = [cos(0.5*nt), ...
        (dth_b/nt)'.*sin(0.5*nt)];
    qP_i2b = quatnormalize(quatmultiply(q_i2b',qdth_b))';
else
    qP_i2b = q_i2b;
end
poseVecP_i2c = poseVec_i2c;
nPoses = size(poseVec_i2c,2);

% Uncertainties
Qa = dt*PSDa*eye(3); % accelerometer noise
Qg = dt*PSDw*eye(3);
Q = blkdiag(Qa,Qg);

% Propagate covariance
% Phi_xx = eye(3,3);
% Phi_xv = dt*eye(3,3);
% Phi_xa = -0.5*(dt^2)*T_i2b_hat'*qf.cpe(am_b);
% Phi_vx = zeros(3);
% Phi_vv = eye(3,3);
% Phi_va = -dt*T_i2b_hat'*qf.cpe(am_b);
% Phi_ax = zeros(3,3);
% Phi_av = zeros(3,3);
% Phi_aa = eye(3,3) - qf.cpe(dth_b);
% Phi = [Phi_xx, Phi_xv, Phi_xa; ...
%     Phi_vx, Phi_vv, Phi_va; ...
%     Phi_ax, Phi_av, Phi_aa];
% if nPts > 0
%     Phi = blkdiag(Phi,eye(3*nPts));
% end
B_xna = -0.5*dt*T_i2b_hat';
B_xng = zeros(3,3);
B_vna = -T_i2b_hat';
B_vng = zeros(3,3);
B_ana = zeros(3,3);
B_ang = eye(3,3);
B = [B_xna, B_xng;...
    B_vna, B_vng; ...
    B_ana, B_ang];

% Augmented covariance matrix
% if nPoses > 0
    P_IMU_IMU = P(1:9,1:9);
    P_IMU_C = P(1:9,10:end); %correlation between errors in IMU state and camera pose estimates
%     P_CC = P(10:end,10:end); %covariance of camera pose estimates
    F_xx = zeros(3); F_xv = eye(3); F_xa = zeros(3); 
    F_vx = zeros(3); F_vv = zeros(3); F_va = -T_i2b_hat'*cpe(am_b);
    F_ax = zeros(3); F_av = zeros(3); F_aa = -cpe(wm_b);
    F = [F_xx, F_xv, F_xa; 
         F_vx, F_vv, F_va; 
         F_ax, F_av, F_aa];
    Phi_IMU_C = expm(F*dt); %xdot = Ax => x = e^(At)x0.
    PP = zeros(size(P));
%     [~,P_hist] = ode45(@(t,P) covOdeFun(t,P,F,B,Q/dt), [0 dt], P_IMU_IMU(:));
%     PP(1:9,1:9) = reshape(P_hist(end,:)',9,9);
    PP(1:9,1:9) = Phi_IMU_C*P_IMU_IMU*Phi_IMU_C' + B*Q*B';
    PP(1:9,10:end) = Phi_IMU_C*P(1:9,10:end); 
    PP(10:end,1:9) = P(10:end,1:9)*Phi_IMU_C'; %(Phi_IMU_C*P_IMU_C)';
    PP(10:end,10:end) = P(10:end,10:end);
% else
%     P_IMU_IMU = P;
%     PP = Phi*P_IMU_IMU*Phi' + B*Q*B';
% end
% L = eig(PP);
% if ~isreal(L) || min(real(L)) < 0
%     L;
% end

end

function [Pdot] = covOdeFun(t,P,F,B,Q)
Pp = reshape(P,9,9);
Pdot = F*Pp + Pp*F' + B*Q*B'; 
Pdot = Pdot(:);
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

function pfHat = lsEstimatePf(pfPkg,poseVec,M)
t_list = pfPkg.t_list;
nt = length(t_list);
pose_inds = [];
meas_inds = [];
for ii = 1:nt
    t = t_list(ii);
    ind = find(M.poseTimestampList == t);
    if ~isempty(ind)
        pose_inds = [pose_inds, ind];
        meas_inds = [meas_inds, ii];
    end
end
pose_list = poseVec(:,pose_inds);
meas_list = pfPkg.meas_list(:,meas_inds);
p1 = pose_list(5:7,1); %first recorded camera position in inertial frame
T_i2c1 = quat2dcm(pose_list(1:4,1)'); %first intertial-to-camera quaternion
p0 = T_i2c1'*[0;0;1] + p1; %guess that the point was 1m in front of the first camera
X0 = p0(1); Y0 = p0(2); Z0 = p0(3);
a0 = X0/Z0; b0 = Y0/Z0; rho0 = 1/Z0;
abr0 = [a0; b0; rho0]; %guess values of alpha, beta, and rho that minimize error in nonlinear measurement
abrHat = lsqnonlin(@(abr)lsFun(abr,meas_list,pose_list(1:4,:), ...
    pose_list(5:7,:)),abr0);
aHat = abrHat(1); bHat = abrHat(2); rhoHat = abrHat(3);
pfHat = (1/rhoHat)*T_i2c1'*[aHat; bHat; 1] + p1
end

function dz = lsFun(abr,z_list,q_list,p_list)
a = abr(1); b = abr(2); rho = abr(3);
n = size(z_list,2); %number of meas-pose pairs for this estimate
p1 = p_list(:,1); %reference everything to first observation
T_i2c1 = quat2dcm(q_list(:,1)'); %quaternions in q_list are q_i2c 
z = z_list(:);
zHat = zeros(size(z));
for ii = 1:n
    p = p_list(:,ii) - p1; %position of camera ii in frame 1
    T = quat2dcm(q_list(:,ii)')*T_i2c1'; %rotation from frame 1 to frame ii
    h = T*[a;b;1] + rho*p;
    zHat(2*ii-1:2*ii) = (1/h(3))*[h(1); h(2)];
end
dz = zHat - z;
end
    



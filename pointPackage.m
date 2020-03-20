function pfPkg = pointPackage(I,tag,meas,t_curr)
% Struct that stores all the data recorded about a single feature point.
% Initialized with first measurement, tag, timestamp, and pose. 
    pfPkg.p_true = I.map(:,tag);
    pfPkg.tag = tag;
    pfPkg.t_list = t_curr; 
    pfPkg.meas_list = meas;
%     pfPkg.pose_list = pose;
end
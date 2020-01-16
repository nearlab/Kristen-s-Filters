function qf = quaternionFunctions()
%QUATERNIONFUNCTIONS Object containing quaternion and general rotation
%operations
qf.cpe = @cpe;
qf.T = @T;
qf.Omega = @Omega;
end

function A = cpe(v)
% Forms a cross-product equivalent matrix from vector v
A = [0, -v(3), v(2);
     v(3), 0, -v(1);
    -v(2), v(1), 0];
end

function T = T(q)
% Makes a passive rotation matrix from quaternion q
qw = q(1); qx = q(2); qy = q(3); qz = q(4);
T = [qw^2+qx^2-qy^2-qz^2, 2*(qx*qy+qw*qz), 2*(qx*qz-qw*qy);
     2*(qx*qy-qw*qz), qw^2-qx^2+qy^2-qz^2, 2*(qy*qz-qw*qx);
     2*(qx*qz-qw*qy), 2*(qy*qz+qw*qx), qw^2-qx^2-qy^2+qz^2];
end

function W = Omega(w)
% Forms Omega matrix for quaternion kinematics from angular velocity w
if size(w,2) ~= 1
    w = w';
end
W = [0, -w';
     w, -cpe(w)];
end



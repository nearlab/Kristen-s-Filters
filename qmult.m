% this function computes the MODIFIED quaternion multiplication
% the right quaternion has the scalar part LAST.
% q3 = qmult(q1,q2)
% the math definition of the quaternion product is
%
% q3(1:3) = q1(4)*q2(1:3) + q2(4)*q1(1:3) - cross(q1(1:3),q2(1:3));
% q3(4)   = q1(4)*q2(4) - dot(q1(1:3),q2(1:3));

function q3 = qmult(q2sf,q1sf)

q1 = [q1sf(2:4); q1sf(1)];  % make it scalar part last
q2 = [q2sf(2:4); q2sf(1)];  % make it scalar part last

q3 = [q1(4)*q2(1) + q2(4)*q1(1) + q1(3)*q2(2) - q1(2)*q2(3);
    q1(4)*q2(2) + q2(4)*q1(2) - q1(3)*q2(1) + q1(1)*q2(3);
    q1(4)*q2(3) + q2(4)*q1(3) - q1(1)*q2(2) + q1(2)*q2(1);
    q1(4)*q2(4) - q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3)];

q3 = [q3(4); q3(1:3)];
end

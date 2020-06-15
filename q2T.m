% function T = q2T(qsf)
% this function computes the 3-by-3 direction cosine matrix T 
% from the right quaternion qsf (scalar part first)
%
% The math formula relating the two is
% T = I - 2*q(4)*S + 2*S^2;
% where
% I = 3x3 identity matrix
% S(v) is the 3x3 cross product skew-symmetric matrix
% S(v) = [0 -v(3) v(2);
%        v(3) 0 -v(1);
%       -v(2) v(1) 0];

function T = q2T(qsf)

r = size(qsf);

if r(1)==1
    qsl = [qsf(2:4)'; qsf(1)];  % make it scalar part last
else
    qsl = [qsf(2:4); qsf(1)];  % make it scalar part last
end

q11 = qsl(1)*qsl(1);
q12 = qsl(1)*qsl(2);
q13 = qsl(1)*qsl(3);
q14 = qsl(1)*qsl(4);
q22 = qsl(2)*qsl(2);
q23 = qsl(2)*qsl(3);
q24 = qsl(2)*qsl(4);
q33 = qsl(3)*qsl(3);
q34 = qsl(3)*qsl(4);
q44 = qsl(4)*qsl(4);

T = [q11-q22-q33+q44 2*(q12+q34) 2*(q13-q24);
     2*(q12-q34) -q11+q22-q33+q44 2*(q23+q14);
     2*(q13+q24) 2*(q23-q14) -q11-q22+q33+q44];


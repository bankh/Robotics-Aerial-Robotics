%% Question 2
% A = rotRes(eye(3), 'z', 0.1)
% B = A * rotRes(A, 'y', 0.1)
% C = B * rotRes(B, 'z', 0.2)
% Symbolic Method
syms ksi the phi
rotZYZ = [cos(ksi) -sin(ksi) 0; sin(ksi) cos(ksi) 0; 0 0 1]*...
         [cos(the) 0 sin(the); 0 1 0; -sin(the) 0 cos(the)]*...
         [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
simplify(rotZYZ)

% Result is
% [ cos(ksi)*cos(phi)*cos(the) - sin(ksi)*sin(phi), - cos(phi)*sin(ksi) - cos(ksi)*cos(the)*sin(phi), cos(ksi)*sin(the)]
% [ cos(ksi)*sin(phi) + cos(phi)*cos(the)*sin(ksi),   cos(ksi)*cos(phi) - cos(the)*sin(ksi)*sin(phi), sin(ksi)*sin(the)]
% [                             -cos(phi)*sin(the),                                sin(phi)*sin(the),          cos(the)]
% the -> cos(the) = 0.995 => acos(0.995) = the = 0.1
%% Question 3
% The coordinates are from Om_b
% a = [a1; a2; a3] -> Sa = [ 0 -a3 a2; a3 0 -a1; -a2 a1 0];
% R_T * R_dot = Om_b -> R_dot = Om_b * inv(R_T)
% R_dot * R_t = Om_s
clear all, clc
R = [ 0.675 -0.1724 0.7174; 0.2474 0.9689 0; -0.6951 0.1775 0.6967]
Om_b = [ 0 -1 0.9689;...
         1  0 -0.2474;...
         -0.9689 0.2474 0];
R_dot = Om_b * inv(R');

Om_s = R_dot * R'

%% Question 4: Axis angle Representations
% u = [a, b, c] and phi = [rad]
% Rot(u,phi) = I cos(phi) + u*u'*(1-cos(phi)) + S(u)*sin(phi)
% tau = r11 + r22 + r33 
% cos(phi) = (tau - 1)/2
% S(u) = (1/(2*sin(phi))) (R - R')
clear all, clc
R = [0.2919 0.643 -0.7081;...
     -0.643 -0.4161 -0.643;...
     -0.7081 0.643 0.2919];
tau = R(1,1) + R(2,2) + R(3,3);
phi = acos((tau-1)/2)
Su = (1/(2*sin(phi)))*(R-R')

%% Question 5: Axis angle Representations
% Quaternion: q = (cos(phi/2), u1*sin(phi/2), u2*sin(phi/2), u3*sin(phi/2))
% Angle of rotation: 2*inv(cos(q0))

% Axis of rotation : u2 = [q1/sqrt(1-q0^2);...
%                          q2/sqrt(1-q0^2);...
%                          q3/sqrt(1-q0^2)]
% where q = [q0, q1, q2, q3]

% 1- Define quaternion: p=(0,p)
% 2- The result after rotation is: p' = q*p*inv(q)=(0,p')
% Compose two rotations q = q2*q1

clear all, clc
R = [-1/3 2/3 -2/3; 2/3 -1/3 -2/3; -2/3 -2/3 -1/3];
quat = qGetQ(R)
u = [quat(2)/sqrt(1-quat(1)^2);...
     quat(3)/sqrt(1-quat(1)^2);...
     quat(4)/sqrt(1-quat(1)^2)]
phi = 2 * acos(quat(1))
% tau = R(1,1) + R(2,2) + R(3,3);
% phi = acos((tau-1)/2)
% Su = (1/(2*sin(phi)))*(R-R')

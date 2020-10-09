clear; clc;
syms q1;
syms q2;
syms q3;

syms d1;
syms a2;
syms d3;

%Forward kinematics symbolical
T = simplify(Rz(q1)*Tz(d1)*Ry(-q2)*Tx(a2)*Tx(q3));

%Forward kinematics test
theta1 = deg2rad(0); theta2 = deg2rad(0); d3 = 10; 
d1 = 10; a2 = 10;

T1 = Rz(theta1)*Tz(d1)*Ry(-theta2)*Tx(a2)*Tx(d3);

%Inverse kinematics
Pos = T1(1:3, 4);
x = Pos(1);
y = Pos(2);
z = Pos(3);

q1_1 = atan2(y,x);
q1_2 = pi + atan2(y,x);
q2_1 = atan2((z-d1), sqrt(x^2 + y^2));
q2_2 = pi - atan2((z-d1), sqrt(x^2 + y^2));
d3_1 = sqrt(x^2 + y^2 + (z-d1)^2) - a2;

%Jacobians classical approach
syms a_2;
syms d_3;
syms d_1;
syms theta_1;
syms theta_2;

T_0 = simplify(Rz(theta_1)*Tz(d_1)*Ry(-theta_2)*Tx(a_2)*Tx(d_3));

px = T_0(1,4);
py = T_0(2,4);
pz = T_0(3,4);

Jv = [diff(px, theta_1), -diff(px, theta_2), diff(px, d_3);
      diff(py, theta_1), -diff(py, theta_2), diff(py, d_3);
      diff(pz, theta_1), -diff(pz, theta_2), diff(pz, d_3); ];
J_om = [0, -sin(theta_1), 0; 0, cos(theta_1), 0; 1, 0, 0];
J_c = simplify([Jv; J_om])
 
%Geometric approach
T0 = eye(4);
T1 = Rz(theta_1) * Tz(d_1);
T2 = T1 * Ry(-theta_2)*Tx(a_2 + d_3);

O0 = T0(1:3, 4);
O1 = T1(1:3, 4);
O2 = T2(1:3, 4);

z0 = T0(1:3, 3);
z1 = T1(1:3, 2);
z2 = T2(1:3, 1);

J1 = [cross(z0,(O2-O0));z0];
J2 = [cross(z1,(O2-O1));z1];
J3 = [z2; 0; 0; 0];

J_G = [simplify(J1), simplify(J2), simplify(J3)]

%Numerical approach

R0 = simplify(T_0(1:3, 1:3));

J1r = Rzd(theta_1)*Tz(d_1)*Ry(-theta_2)*Tx(a_2)*Tx(d_3)*...
    [R0^-1 zeros(3,1); 0 0 0 1];
J1 = [J1r(1,4) J1r(2,4) J1r(3,4) J1r(3,2) J1r(1,3) J1r(2,1)]';

J3r = Rz(theta_1)*Tz(d_1)*Ryd(-theta_2)*Tx(a_2)*Tx(d_3)*...
    [R0^-1 zeros(3,1); 0 0 0 1];
J3 = [J3r(1,4) J3r(2,4) J3r(3,4) J3r(3,2) J3r(1,3) J3r(2,1)]';

J5r = Rz(theta_1)*Tz(d_1)*Ry(-theta_2)*Tx(a_2)*Txd(d_3)*...
    [R0^-1 zeros(3,1); 0 0 0 1];
J5 = [J5r(1,4) J5r(2,4) J5r(3,4) J5r(3,2) J5r(1,3) J5r(2,1)]';

J_N = [simplify(J1), simplify(J3), simplify(J5)]

%Singularities

det_Jv = simplify(det(Jv))

%Velocity
syms Theta1 Theta2 D3 t 

Theta1 = sin(t);
Theta2 = cos(2*t);
D3 = sin(3*t);

J_time = simplify(subs(J_G, {d_1, a_2, theta_1, theta_2, d_3}, {10, 10, Theta1, Theta2, D3}));
q = [Theta1 Theta2 D3]';
qdot = diff(q);

xi = simplify(J_time * qdot);

time = 0:0.01:10;
xi_time = subs(xi, {t}, {time});
xi_t = double(xi_time);

figure(1);
plot(time, xi_t(1:3, :), 'LineWidth', 2)
xlabel('Time of moving')
ylabel('Linear velocity value')
title('Linear velocities plots')
legend('v_x', 'v_y', 'v_z', 'Fontsize', 15)
grid on

figure(2);
plot(time, xi_t(4:6, :), 'LineWidth', 2)
xlabel('Time of moving')
ylabel('Angular velocity value')
title('Angular velocities plots')
grid on
legend('\omega_x', '\omega_y', '\omega_z', 'Fontsize', 15)



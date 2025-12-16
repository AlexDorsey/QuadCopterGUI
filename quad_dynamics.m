function [xdot, accel_body_for_filter, angularVelocity_body] = quad_dynamics(x, moment, force, params)
    %x(1): phi: pitch
    %x(2): theta: roll
    %x(3): psi: yaw
    %x(4:6): body rates
    %x(7): east
    %x(8): north
    %x(9): up
    %x(10): v_east
    %x(11): v_north
    %x(12): v_up
    J = params.J;
    m = params.mass_kg;
    
    xdot = x*0;
    Eul_Angles = x(1:3,1);
    phi     = Eul_Angles(1);
    theta   = Eul_Angles(2);
    psi     = Eul_Angles(3);
    omega_QE = x(4:6,1);
    
    xdot(1:3,1) = S_Phi_Theta_inv(phi,theta)*omega_QE;
    xdot(4:6,1) = J\(-cross(omega_QE,J*omega_QE)+moment);
    
    xdot(7:9,1) = x(10:12,1);
    R_enu_to_body = get_R_enu_to_body(phi, theta, psi);
    R_body_to_enu = transpose(R_enu_to_body);

    V_dot_enu = [0;0;-9.81] + R_body_to_enu*[0;0;force/m];
    % [O3psi, O2theta, O1phi] = orientation_matrices(phi, theta, psi);
    % V_dot_enu = [0;0;-9.81] + force/m*transpose(O3psi) * transpose(O2theta) * transpose(O1phi)*[0;0;1];
    xdot(10:12,1) = V_dot_enu;
     
    accel_body_for_filter = R_enu_to_body * [0;0;9.81+force/m];
    angularVelocity_body = omega_QE;

end


function R_enu_to_body = get_R_enu_to_body(phi, theta, psi)
    % Rotation from ENU (world) to body using ZYX (yaw-psi, pitch-theta, roll-phi).
    % v_body = R * v_enu
    
    cphi = cos(phi);   sphi = sin(phi);
    cth  = cos(theta); sth  = sin(theta);
    cps  = cos(psi);   sps  = sin(psi);
    
    Rz = [cps -sps 0;
          sps  cps 0;
          0     0  1];
    
    Ry = [ cth 0 sth;
            0  1  0;
          -sth 0 cth];
    
    Rx = [1   0    0;
          0 cphi -sphi;
          0 sphi  cphi];
    
    R_enu_to_body = (Rz*Ry*Rx)';   % transpose of body->world
end

function S = S_Phi_Theta_inv(Phi,Theta)
S = [1 sin(Phi)*tan(Theta) cos(Phi)*tan(Theta);
    0 cos(Phi) -sin(Phi);
    0 sin(Phi)*sec(Theta) cos(Phi)*sec(Theta)];
end

function [O3psi, O2theta, O1phi] = orientation_matrices(phi, theta, psi)
    O3psi = [cos(psi) sin(psi) 0;
        -sin(psi) cos(psi) 0;
        0 0 1];
    O2theta = [cos(theta) 0 -sin(theta);
        0 1 0
        sin(theta) 0 cos(theta)];
    O1phi = [1 0 0
        0 cos(phi) sin(phi);
        0 -sin(phi) cos(phi) ;];
end



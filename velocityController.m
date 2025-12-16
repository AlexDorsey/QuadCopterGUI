function [force_out, angleCommands] = velocityController(iteration, x, velocityCommands, dt, psi_cmd, USER_PARAMS)

    persistent error_v_e_int error_v_n_int error_v_u_int phi_cmd_filt theta_cmd_filt
    if iteration == 1 || isempty(error_v_e_int)
        error_v_e_int = 0;
        error_v_n_int = 0; 
        error_v_u_int = 0;
        phi_cmd_filt = 0;
        theta_cmd_filt = 0; 
    end
    % extract states
    phi  = x(1); 
    theta= x(2); 
    psi  = x(3);
    v_e  = x(10);
    v_n  = x(11);
    v_u  = x(12);
    
    % roation matricies:
    R_enu_to_body = get_R_enu_to_body(phi, theta, psi);
    % commands
    v_e_cmd = velocityCommands(1);
    v_n_cmd = velocityCommands(2);
    v_u_cmd = velocityCommands(3);
    
    % errors
    VelocityErrorLimit = 2.5;
    error_v_e = clip(v_e_cmd - v_e, -VelocityErrorLimit, VelocityErrorLimit);
    error_v_n = clip(v_n_cmd - v_n, -VelocityErrorLimit, VelocityErrorLimit);
    error_v_u = clip(v_u_cmd - v_u, -VelocityErrorLimit, VelocityErrorLimit);
    
    % error ints
    error_v_e_int = error_v_e_int + error_v_e*dt;
    error_v_n_int = error_v_n_int + error_v_n*dt;
    error_v_u_int = error_v_u_int + error_v_u*dt;

    % constants
    g = 9.81;          % m/s^2
    m = USER_PARAMS.mass_kg;           % kg 
    
    % desired force:
    Kp = 0.75*m; Ki = 0.01*m;
    force_east = Kp*error_v_e + Ki*error_v_e_int;
    force_north = Kp*error_v_n + Ki*error_v_n_int;
    force_up = Kp*error_v_u + Ki*error_v_u_int + m*g;
    
    % desired forces vectors
    force_enu = [force_east; force_north; force_up];
    force_out = norm(force_enu);

    if force_out < 1 % clipping at 1 Newton
        force_out = 1;
    end
    
    [phi_cmd, theta_cmd] = force_to_roll_pitch_yaw_corrected(force_enu, psi_cmd);

    % clip roll and pitch, wrap yaw
    lim = 65*pi/180;
    phi_cmd   = clip(phi_cmd,   -lim, lim);
    theta_cmd = clip(theta_cmd, -lim, lim);
    
    theta_cmd_filt = firstOrderFilter(theta_cmd_filt, theta_cmd, dt,  USER_PARAMS.AngleCommandFilterTau);
    phi_cmd_filt = firstOrderFilter(phi_cmd_filt, phi_cmd, dt,  USER_PARAMS.AngleCommandFilterTau);
    angleCommands = [phi_cmd_filt; theta_cmd_filt; wrapToPi(psi_cmd)];

    % for I = 1:10000
    %     command = 10*sign(sin(I*dt));
    %     if I ==1 
    %         response = command;
    %     end
    %     [response] = firstOrderFilter(response, command, dt, 0.5);
    %     responseList(I) = response;
    %     commandList(I) = command;
    % end
    % 
    % figure('Position',[10,10,1000,1000]); plot(commandList,'k--'); hold on; plot(responseList); drawnow;
end


function [phi_cmd, theta_cmd] = force_to_roll_pitch_yaw_corrected(force_enu, psi_cmd)

    % normalize
    fn = norm(force_enu);
    if fn < 1e-12
        phi_cmd = 0;
        theta_cmd = 0;
        return
    end
    fhat = force_enu / fn;

    % rotate world force by -yaw so extraction assumes yaw = 0
    c = cos(-psi_cmd); s = sin(-psi_cmd);
    Rz_mpsi = [c -s 0;
               s  c 0;
               0  0 1];
    f0 = Rz_mpsi * fhat;    

    phi_cmd   = -asin( clip(f0(2), -1, 1) );
    theta_cmd = atan2( f0(1), f0(3) );
end

function R_enu_to_body = get_R_enu_to_body(phi, theta, psi)
    
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
    
    R_enu_to_body = (Rz*Ry*Rx)';   
end

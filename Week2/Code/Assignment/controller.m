function [ u ] = pd_controller(~, s, s_des, params)
%PD_CONTROLLER  PD controller for the height
%
%   s: 2x1 vector containing the current state [z; v_z]
%   s_des: 2x1 vector containing desired state [z; v_z]
%   params: robot parameters

u = 0;

K_p = 200; K_v = 20;
e = s_des - s;
u = params.mass*(K_p*e(1)+K_v*e(2)+params.gravity);
if(u>params.u_max) u = params.u_max; end
if(u<params.u_min) u = params.u_min; end

end


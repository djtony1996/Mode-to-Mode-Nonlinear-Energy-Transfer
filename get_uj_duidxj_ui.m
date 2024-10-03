% This code is to calculate the term:
%               u_j * \frac{\partial u_i}{\partial x_j} * u_i
function [nineterms] = get_uj_duidxj_ui(u_F1,v_F1,w_F1,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,u_F2,v_F2,w_F2)

nineterms    = u_F1 .* duF_dx .* u_F2 + u_F1 .* dvF_dx .* v_F2 + u_F1 .* dwF_dx .* w_F2 ...
             + v_F1 .* duF_dy .* u_F2 + v_F1 .* dvF_dy .* v_F2 + v_F1 .* dwF_dy .* w_F2 ...
             + w_F1 .* duF_dz .* u_F2 + w_F1 .* dvF_dz .* v_F2 + w_F1 .* dwF_dz .* w_F2;

end



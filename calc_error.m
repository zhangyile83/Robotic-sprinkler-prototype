function [ error ] = calc_error( target_x, target_y, p1, p2, v0, gama, alpha, beta, diameter, wind_speed)

[drop_x,drop_y,~]=droppoint(v0*p1, gama*p2, alpha, beta, diameter, wind_speed);
error= ((target_x-drop_x)^2+(target_y-drop_y)^2)^0.5;

end
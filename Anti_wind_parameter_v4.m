%%
%__ Multi-time corrector for both velocity and angle
%__ NON-Debug Mode

% _____________ Input parameters ________________
% v0 is the flow velocity without wind
% Dx is the displacement at x-coordinate
% Dy is the displacement at y-coordinate
% w is the wind speed
% gamaw is the horizontal angle of the wind
% MeanD is the mean diameter
% t is the dropping time of case without correction

% _____________ Output parameters ________________
% alpha1 is the velocity correction parameter
% alpha2 is the angle correction parameter
% t1 is the dropping time of the case with optimal correction

% ____________find the dropping time for correction case___________________
% find the time that the droplet touch the ground after correction

%%
function [output_p1, output_p2] = Anti_wind_parameter_v4(v0, target_x, target_y, gama, diameter, wind_speed, beta)

%%
%constants
increase=1;%decrease=-1;
threshold = 0.05;
scale_p1 = 0.05;
scale_p2 = 0.002;
alpha = 30;

persistent p1;
persistent p2;
persistent adjust_direction_p1;
persistent adjust_direction_p2;
if isempty(p1)
    p1=1;
end
if isempty(p2)
    p2=1;
end
if isempty(adjust_direction_p1)
    adjust_direction_p1=increase;
end
if isempty(adjust_direction_p2)
    adjust_direction_p2=increase;
end
error_min=calc_error(target_x, target_y, p1, p2, v0, gama, alpha, beta, diameter, wind_speed);

%%
while(1)
    %adjust p1
    better=0;
    error_tmp=calc_error(target_x, target_y, p1+adjust_direction_p1*scale_p1, p2, v0, gama, alpha, beta, diameter, wind_speed);
    if error_tmp<error_min%如果这个方向，则一直照这个方向运行直到极点
        better=1;
    else%如果这个方向不能优化，则尝试另一个方向
        error_tmp=calc_error(target_x, target_y, p1-adjust_direction_p1*scale_p1, p2, v0, gama, alpha, beta, diameter, wind_speed);
        if error_tmp<error_min
            adjust_direction_p1=-adjust_direction_p1;
            better=1;
        end
    end
    while better %由于matlab没有dowhile，所以这么写
        p1=p1+adjust_direction_p1*scale_p1;
        error_min=error_tmp;
        error_tmp=calc_error(target_x, target_y, p1+adjust_direction_p1*scale_p1, p2, v0, gama, alpha, beta, diameter, wind_speed);
        if error_tmp>error_min
            break;
        end
    end
    
    %%
    %adjust p2
    better=0;
    error_tmp=calc_error(target_x, target_y, p1, p2+adjust_direction_p2*scale_p2, v0, gama, alpha, beta, diameter, wind_speed);
    if error_tmp<error_min%如果这个方向，则一直照这个方向运行直到极点
        better=1;
    else%如果这个方向不能优化，则尝试另一个方向
        error_tmp=calc_error(target_x, target_y, p1, p2-adjust_direction_p2*scale_p2, v0, gama, alpha, beta, diameter, wind_speed);
        if error_tmp<error_min
            adjust_direction_p2=-adjust_direction_p2;
            better=1;
        end
    end
    while better %由于matlab没有dowhile，所以这么写
        p2=p2+adjust_direction_p2*scale_p2;
        error_min=error_tmp;
        error_tmp=calc_error(target_x, target_y, p1, p2+adjust_direction_p2*scale_p2, v0, gama, alpha, beta, diameter, wind_speed);
        if error_tmp>error_min
            break;
        end
    end
    if error_min>threshold
        scale_p1 = scale_p1/2;
        scale_p2 = scale_p2/2;
    else
        break;
    end
end
%%
%debugging infomation
persistent error_max;
persistent error_total;
persistent number;
if isempty(error_max)
    error_max=0;
end
if isempty(error_total)
    error_total=0;
end
if isempty(number)
    number=0;
end
number=number+1;
error_total=error_total+error_min;
if(error_min>error_max)
    error_max=error_min;
end
% fprintf('average error: %f\n', error_total/number);
% fprintf('max     error: %f\n', error_max);
% if error_min>threshold
%     fprintf('Exceeding threshold!');
%     pause(10);
%     return;
% end

%%
%assign to output
output_p1=p1;
output_p2=p2;
end

%%%%%%%%%%%%%%%%%%%%%% Wind shifting Technique %%%%%%%%%%%%%%%%%%%%%%%%
% __ wind effect correction for 360 degree case
% __ Multi droplet diameter
% __ Multiple spraying circle
% __ corrections compute in advance
% __ Smooth version
%clear;
%%
%input
alpha = 30;  % The angle of elevation of the nozzle
MeanD = 0.01;  % Droplet mean diameter(should be the value defined in the range)
vertex_x=[23;-23;-23;23];
vertex_y=[23;-23;23;-23];
distance_proportion_start = 0.5; 
distance_proportion_step_length = 0.05;
distance_proportion_end = 1;
wind_speed_start = 0.447*3;
wind_speed_step_length = 0.447;
wind_speed_end = 6.7050;
beta_start = 15; 
beta_step_length = 5;
beta_end = 360;
angle_step_length=10;
%%
%debugging information
count=0;
TotalTime=0;
TotalParaTime=0;
%%
%initialization
A1 = zeros(16, 72, 181, 23);%16个风速，72个风向，181转角，23米半径
A2 = zeros(16, 72, 181, 23);
%%
%calcylate contour information
[contour_distance, angle_start, angle_end]=calc_contour(vertex_x, vertex_y, angle_step_length);
index_circle=1;
for distance_proportion=distance_proportion_start:distance_proportion_step_length:distance_proportion_end        
    index_wind_speed = 1;
    for wind_speed = wind_speed_start:wind_speed_step_length:wind_speed_end%风速，调试用for wind_speed = 0:0.447:6.7050;        
        index_beta = 1;
        for beta = beta_start:beta_step_length:beta_end;%风向，调试用for beta = 0:5:359;
            index_gama = 1; % 181 in total
            %_______ Spraying loop at different angle in each spraying circle ______
            for gama=angle_start:angle_step_length:angle_end
                ticMark_total=tic;%debugging information
                distance=contour_distance(index_gama)*distance_proportion;
                v0=Characteristic(alpha, distance, MeanD); %Find the required flow speed for different spraying angles according
                ticMark_para=tic;
                [p1, p2] = Anti_wind_parameter_v4(v0, distance*cosd(gama), distance*sind(gama), gama, MeanD, wind_speed, beta);%_____________Find shifting parameters_________________
                A1(index_wind_speed, index_beta, index_gama, index_circle) = p1;
                A2(index_wind_speed, index_beta, index_gama, index_circle) = p2;
                %%                
                %deubgging information                
                TotalParaTime=TotalParaTime+toc(ticMark_para);
                index_gama = index_gama + 1;
                count = count + 1;
                TotalTime=TotalTime+toc(ticMark_total);
%                 display(count);
%                 fprintf('param: %f\n',TotalParaTime/count);
%                 fprintf('total: %f\n',TotalTime/count);
                fprintf('%d: %f %f\n',count, p1*v0, gama*p2-360);                
                control_series(count,1)=p1*v0;
                control_series(count,2)=gama*p2;                
            end
            index_beta = index_beta + 1;
        end
        index_wind_speed = index_wind_speed + 1;
    end    
    index_circle=index_circle+1;
end
function [ contour_distance, angle_start, angle_end ] = calc_contour( vertex_x, vertex_y, angle_scale )
vertex_number=size(vertex_x,1);
%%
%calculate the angle
angle=zeros(vertex_number,1);
for i=1:vertex_number
    [angle(i,1),~]=cart2pol(vertex_x(i),vertex_y(i));
    if(angle(i,1)<0)
        angle(i,1)=angle(i,1)+2*pi;
    end
    angle(i,1)=round(angle(i,1)*180/pi);
end
vertexes = sortrows([vertex_x,vertex_y,angle],3);
%%
%judge if it's on a vertex
on_a_vertex=0;
for i=1:vertex_number
    if vertex_x(i)==0 && vertex_y(i)==0
        on_a_vertex =i;        
        on_the_line=1;
        vertex_number=vertex_number-1;
        vertexes(i,:)=[];
        break;
    end
end
%%
%judge if it's on the line
if(on_a_vertex==0)
    on_the_line=0;
    for i=1:vertex_number-1
        if rank([[0,0];vertexes(i,1:2);vertexes(i+1,1:2)])==1
            on_the_line=i;
            break;
        end
    end
    if rank([[0,0];vertexes(vertex_number,1:2);vertexes(1,1:2)])==1
        on_the_line=vertex_number;
    end
end
%%
%assign angle
if on_the_line==0
    angle_start=0;
    angle_end=0:angle_scale:359;
    angle_end=angle_end(end);
    
else
    angle_start=min(angle);
    angle_end=0:angle_scale:max(angle);
    angle_end=angle_end(end);    
end
%%
%calculate contour_distance
count_tmp=1;

if on_the_line==0
    contour_distance=zeros(1,size(0:angle_scale:359,2)); 
    angle_tmp=0;
    while (angle_tmp<vertexes(1,3))
        contour_distance(count_tmp)=calc_intersection(angle_tmp,vertexes(1,1),vertexes(1,2),vertexes(vertex_number,1),vertexes(vertex_number,2));
        count_tmp=count_tmp+1;
        angle_tmp=angle_tmp+angle_scale;
    end
    for index_vertex=1:vertex_number-1
        while (angle_tmp<vertexes(index_vertex+1,3))        
            contour_distance(count_tmp)=calc_intersection(angle_tmp,vertexes(index_vertex,1),vertexes(index_vertex,2),vertexes(index_vertex+1,1),vertexes(index_vertex+1,2));
            count_tmp=count_tmp+1;
            angle_tmp=angle_tmp+angle_scale;
        end
    end
    while (angle_tmp<360)        
        contour_distance(count_tmp)=calc_intersection(angle_tmp,vertexes(1,1),vertexes(1,2),vertexes(vertex_number,1),vertexes(vertex_number,2));
        count_tmp=count_tmp+1;
        angle_tmp=angle_tmp+angle_scale;
    end
else     
    contour_distance=zeros(1,size((angle_start:angle_scale:angle_end),2));
    angle_tmp=angle_start;
    for index_vertex=1:vertex_number-1        
        while (angle_tmp<vertexes(index_vertex+1,3))        
            contour_distance(count_tmp)=calc_intersection(angle_tmp,vertexes(index_vertex,1),vertexes(index_vertex,2),vertexes(index_vertex+1,1),vertexes(index_vertex+1,2));
            count_tmp=count_tmp+1;
            angle_tmp=angle_tmp+angle_scale;
        end
    end
end
angle_start=angle_start+360;
angle_end=angle_end+360;
end
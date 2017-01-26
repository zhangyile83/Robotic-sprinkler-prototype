function [distance]= calc_intersection( angle,x1,y1,x2,y2 )
point1_1=[0,0];
point1_2=[1*cosd(angle),sind(angle)];
point2_1=[x1,y1];
point2_2=[x2,y2];
if point1_1(1)==point1_2(1)
    X=point1_1(1);
    k2=(point2_2(2)-point2_1(2))/(point2_2(1)-point2_1(1));
    b2=point2_1(2)-k2*point2_1(1);
    Y=k2*X+b2;
end
if point2_1(1)==point2_2(1)
    X=point2_1(1);
    k1=(point1_2(2)-point1_1(2))/(point1_2(1)-point1_1(1));
    b1=point1_1(2)-k1*point1_1(1);
    Y=k1*X+b1;
end
if point1_1(1)~=point1_2(1)&point2_1(1)~=point2_2(1)
    k1=(point1_2(2)-point1_1(2))/(point1_2(1)-point1_1(1));
    k2=(point2_2(2)-point2_1(2))/(point2_2(1)-point2_1(1));
    b1=point1_1(2)-k1*point1_1(1);
    b2=point2_1(2)-k2*point2_1(1);
    if k1==k2
        X=[];
        Y=[];
    else
        X=(b2-b1)/(k1-k2);
        Y=k1*X+b1;
    end
end
distance=sqrt(X^2+Y^2);
end
%%Here is an example of calculating Rotate Matrix and Position Matrix of
%%Camera
%Projection Matrix
P=[0.559445 -0.775775 -0.291879 -0.831389;
-0.722298 -0.283554 -0.630780 -0.454073;
0.406580 0.563711 -0.718974 1.121994];

%Original rotate matrix
oriR=[0.559445 -0.775775 -0.291879;
-0.722298 -0.283554 -0.630780;
0.406580 0.563711 -0.718974];

%calculate the inverse matrix of oriR
R=inv(oriR)
T1=[-0.831389;-0.454073;1.121994];
position=-1.0*R*T1



%calculate Euler Angle according to Rotate Matrix
yd = asind(-R(3,1));
if cos(yd) ~= 0
    if R(1,1)/cosd(yd) > 1
        zd = acosd(1);
    elseif R(1,1)/cosd(yd) < -1
        zd = acosd(-1);
    else
        zd = atand(R(2,1)/R(1,1));
    end

    if R(3,3)/cosd(yd) > 1
        xd = acosd(1);
    elseif R(3,3)/cosd(yd) < -1
        xd = acosd(-1);
    else
        xd = asind(R(3,2)/cosd(yd));
    end
else

end

%output: [-pi/2,pi/2]
r = [xd; yd; zd]

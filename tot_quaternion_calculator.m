function quaternion = tot_quaternion_calculator(t)
lend3 = sqrt(25 + 4*cos(3*t)^2+16*cos(3*t));
R11 = cos(3*t)*cos(2*t);
R12 = (4*cos(2*t)*sin(3*t)+2*cos(3*t)*cos(2*t)*sin(3*t)-3*sin(2*t)*sin(3*t)^2-3*cos(3*t)*sin(3*t)) / lend3;
R13 = (-4*sin(2*t)-3*sin(3*t)*cos(2*t)-2*cos(3*t)*sin(2*t))/ lend3;
R21 = cos(3*t)*sin(2*t);
R22 = ((3*cos(3*t)^2*cos(2*t)+2*sin(2*t)*sin(3*t)+cos(3*t)*sin(2*t)*sin(3*t)+3*cos(2*t)*sin(3*t)^2)) / lend3;
R23 = (4*cos(2*t)-3*sin(3*t)*sin(2*t)+2*cos(3*t)*cos(2*t)) / lend3;
R31 = sin(3*t);
R32 = ((-2*sin(2*t)^2*cos(3*t)-cos(3*t)^2*sin(2*t)^2-4*cos(3*t)*cos(2*t)^2-2*cos(3*t)^2*cos(2*t)^2)) / lend3;
R33 = 3*cos(3*t)/lend3;
qw = ((1 + R11 + R22 + R33)^(1/2))/2;
if t > pi
    qw = -qw;
end
qx = (R32 - R23)/(4*qw);
qy = (R13 - R31)/(4*qw);
qzed = (R21 - R12)/(4*qw);
lenquaternion = (qw^2+qx^2+qy^2+qzed^2)^(1/2);
quaternion = [qx/lenquaternion,qy/lenquaternion,qzed/lenquaternion,qw/lenquaternion];
end
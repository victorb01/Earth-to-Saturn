function [tim,dtimdz,capA,capY,iflag] = getTimeLambert(z,mu,k,theta,r0,r,c,s,dc,ds)
%UNTITLED Summary of this function goes here
%   tim is the time of flight as a function of z
%   dtimdz is derivative of time wrt z
%   capA,capY are intermediate variables
%   iflag is 0 on normal return, if -1 then Y is negative and so the z is
%   in an invalid region (too negative)
%
%   z is the guess variable
%   {mu,k,theta,r0,r} are parameters of the problem 
%   c,s,dc,ds are the C(z) and S(z) and their first and second derivatives
%   
tim=0;
dtimdz=0;

iflag=0;

t1 = cos(theta);
t5 = sqrt(r0 * r * (0.1e1 + t1));
t6 = t5 * k;
capA = t6;

t8 = -s * z + 1;
t9 = sqrt(c);
t10 = 0.1e1 / t9;
capY = -t10 * t8 * t6 + r + r0;
if(capY<0)
    iflag=-1;
    return
end


t13 = sqrt(mu);
t14 = 0.1e1 / t13;
t15 = 0.1e1 / c;
t16 = t15 * capY;
t17 = sqrt(t16);
t18 = t17 * t16;
t20 = sqrt(capY);
tim = (s * t18 + t20 * t6) * t14;


t23 = s * t17;
t25 = -ds * z - s;
t29 = 0.1e1 / t9 / c;
t30 = t29 * t8;
t34 = -t10 * t25 * capA + dc * t30 * capA / 0.2e1;
t36 = c ^ 2;
t37 = 0.1e1 / t36;
t38 = t37 * capY;
t40 = -dc * t38 + t15 * t34;
t44 = 0.1e1 / t20;
dtimdz = (0.3e1 / 0.2e1 * t40 * t23 + ds * t18 + t34 * t44 * capA / 0.2e1) * t14;




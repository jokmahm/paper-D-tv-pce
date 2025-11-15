function y = Jump(a,z1,z2,mod)

switch mod
    case 1
% Example 1
if z2<=0.5 
    y = mortality(a); 
else
    y = 0; 
end

    case 2
% Example 2
if z2<=0.5+0.1*sin(2*pi*z1)
   y = mortality(a); 
else
   y = 0;
end

    case 3
% Example 3
r1 = 0.8;
r2 = 0.8;

if (z1-r1)^2+(z2-r2)^2 <= 1/10
    y = 0;
else
    y = mortality(a);
end
end
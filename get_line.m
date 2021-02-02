function [] = get_line(x1,y1,x2,y2)

m = (y1-y2)/(x1-x2);
b = -m*x1 +y1;
fprintf("y = %0.3gx + %0.3g\n", m, b)
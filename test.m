clc
clear all
addpath('fun')

y = 1:100;
y = y';
x = makeYlag(y, 3);

G = [cos(pi/3), -sin(pi/3);
     sin(pi/3), cos(pi/3)];
matexp(G, 0)
matexp(G, 1) 
matexp(G, 3)
matexp(G, 6)
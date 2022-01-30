clc
clear all
addpath('fun')

y = 1:100;
y = y';
x = makeYlag(y, 3);
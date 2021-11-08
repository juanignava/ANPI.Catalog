clc;clear;
pkg load symbolic
syms x;
a=2; b=5;
f='log(x)';
fs=sym(f);
y=((b-a)*x+(b+a))/2;
gs=(b-a)/2*subs(fs,x,y)
gn=matlabFunction(gs)
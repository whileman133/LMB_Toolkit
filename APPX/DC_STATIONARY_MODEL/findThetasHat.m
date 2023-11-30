% findThetasHat.m

clear; close all; clc;
syms s Ds dUocp Cs tauf nf;

zeta = (s/Ds)*((1+s*tauf)/tauf/Ds)^(nf-1);
beta = sqrt(zeta);
h = (beta^2+3*(1-beta*coth(beta)))/beta^2/(1-beta*coth(beta));
h0 = limit(h,s,0);
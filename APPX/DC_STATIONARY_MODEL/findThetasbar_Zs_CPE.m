% findThetasbar.m

clear; close all; clc;
syms s x;
syms R0 R2 R3 V0 Cs Rdl Cdl Rf tauf nf Ds ThetasBar;
syms sp kp W;

% Definition of R3 to be substituted later.
dfnR3 = ( 2/sp+3/kp-(1/sp+1/kp)*(x+2)/2-W*(x-4)/2/kp )*(x-2);

% Compute dc integrator-removed solid-liquid resistance.
Zk = R0 + V0*ThetasBar + 1/Cs/(1+s*tauf)^(nf-1)/s;
Zdl = Rdl + 1/Cdl/s;
Zse = Zk*Zdl/(Zk+Zdl)+Rf;
res0 = limit(s*Zse,s,0);
res0 = simplify(res0);
ZseStar = Zse - res0/s;
ZseStar = collect(simplify(expand(ZseStar)),s);
RseStar = limit(ZseStar,s,0);

% Solve for ThetasBar and R2.
solnThetasBar = solve(R2+R3==-RseStar,ThetasBar);
solnThetasBar = simplify(solnThetasBar);
solnThetasBarSubR3 = subs(solnThetasBar,R3,dfnR3);
solnR2 = solve(int(solnThetasBarSubR3,x,2,3)==0,R2);
solnR2 = simplify(solnR2);
% findThetasbar.m
% Includes both double-layer and solid-diffusion CPEs!

clear; close all; clc;
syms s x;
syms R0 R2 R3 V0 Cs Rdl Cdl Rf tauf nf Ds taudl ndl ThetasBar;
syms sp kp W dUocp Rct;
syms Cdldc Dsdc;

% Definitions to be substituted later.
dfnV0 = -dUocp*(1+Cdldc/Cs);
dfnV0ALT = -dUocp*(1+Cdldc/Cs)^2;
dfnR0 = Rct + 1/15/Dsdc/Cs;
dfnR3 = ( 2/sp+3/kp-(1/sp+1/kp)*(x+2)/2-W*(x-4)/2/kp )*(x-2);
dfnThetasHat = -1/15/Dsdc/dUocp/(Cs+Cdldc);

% Compute dc integrator-removed solid-liquid resistance.
Zk = R0 + V0*ThetasBar + (1+s*tauf)^(1-nf)/Cs/s;
Zdl = Rdl + ((Cdl/taudl)*(1+taudl*s))^(1-ndl)/Cdl/s;
Zse = Zk*Zdl/(Zk+Zdl)+Rf;
res0 = limit(s*Zse,s,0);
res0 = simplify(res0);
ZseStar = Zse - res0/s;
ZseStar = collect(simplify(expand(ZseStar)),s);
RseStar = limit(ZseStar,s,0);
RseStar = simplify(RseStar,'Steps',1000);

% Solve for ThetasBar and R2.
solnThetasBar = solve(R2+R3==-RseStar,ThetasBar);
solnThetasBar = simplify(solnThetasBar,'Steps',1000);
solnThetasBarSubR3 = subs(solnThetasBar,R3,dfnR3);
solnR2 = solve(int(solnThetasBarSubR3,x,2,3)==0,R2);
solnR2 = simplify(solnR2,'Steps',1000);

% Substutions and simplification.
solnThetasBar1 = subs(solnThetasBar,taudl*(Cdl/taudl)^ndl,Cdldc);
solnThetasBar1 = subs(solnThetasBar1,V0,dfnV0);
solnThetasBar1 = subs(solnThetasBar1,R0,dfnR0);
solnThetas1 = dfnThetasHat + solnThetasBar1;
solnThetas1 = simplify(solnThetas1,'Steps',1000);
solnThetas1 = subs(solnThetas1,taudl*(Cdl/taudl)^ndl,Cdldc);
solnThetas1 = simplify(solnThetas1,'Steps',1000);
solnThetasBar2 = subs(solnThetasBar,V0,dfnV0ALT);
solnThetasBar2 = subs(solnThetasBar2,R0,dfnR0);
solnThetas2 = dfnThetasHat + solnThetasBar2;
solnThetas2 = simplify(solnThetas2,'Steps',1000);
solnThetas2 = subs(solnThetas2,taudl*(Cdl/taudl)^ndl,Cdldc);
solnThetas2 = simplify(solnThetas2,'Steps',1000);

solnR2_1 = subs(solnR2,taudl*(Cdl/taudl)^ndl,Cdldc);
solnR2_1 = subs(solnR2_1,taudl^2*(Cdl/taudl)^(2*ndl),Cdldc^2);
solnR2_1 = subs(solnR2_1,(Cdl/taudl)^ndl,Cdldc/taudl);
solnR2_1 = simplify(solnR2_1,'Steps',1000);
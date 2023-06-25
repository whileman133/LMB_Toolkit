function x = denormalizePosition(xtilde,p2dm)
%DENORMALIZEPOSITION Convert dimensionless position vector back to meters.

x = zeros(size(xtilde));
dll = 0<=xtilde&xtilde<1;
sep = 1<=xtilde&xtilde<2;
pos = 2<=xtilde&xtilde<=3;
p = getCellParams(p2dm,'*.L');
x(dll) = xtilde(dll)*p.Ldll;
x(sep) = p.Ldll + (xtilde(sep)-1)*p.Lsep;
x(pos) = p.Ldll + p.Lsep + (xtilde(pos)-2)*p.Lpos;

end
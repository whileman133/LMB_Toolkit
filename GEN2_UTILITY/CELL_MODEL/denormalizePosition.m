function x = denormalizePosition(xtilde,p2dm)
%DENORMALIZEPOSITION Convert dimensionless position vector back to meters.

x = zeros(size(xtilde));
dll = 0<=xtilde&xtilde<1;
sep = 1<=xtilde&xtilde<2;
pos = 2<=xtilde&xtilde<=3;
p = getCellParam(p2dm,'dll.L sep.L pos.L');
x(dll) = xtilde(dll)*p.dll.L;
x(sep) = p.dll.L + (xtilde(sep)-1)*p.sep.L;
x(pos) = p.dll.L + p.sep.L + (xtilde(pos)-2)*p.pos.L;

end
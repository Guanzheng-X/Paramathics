
% declare all the symbols which you want to track as symbol
a = sym('A');
qf = sym('Qf');
pi = sym('Pi');
pit = sym('PiT');
ri = sym('Ri');
xi = sym('Xi');
b = sym('B');

% declare q transpose so avoid cancellation
qft = sym('QfT');
qlt = sym('QlT');
qct = sym('QcT');
rlt = sym('RlT');
rlt1 = sym('RlT1');
pft = sym('PfT');
plt = sym('PlT');

% while defining alp make sure that you use different q then qi
% because it will cancel out, while we update x and r
% this code uses qit -> "QiT"


% iteration of cg when bit is flipped in a and qf gets corrupted
%qf = a * pi;
alpf = (ri^2) / (pit *qft);
xf = xi + (alpf * pi);
rf = ri - (alpf * qf);
btf = (rf^2) / (ri^2);
pf = rf + (btf * pi);

% iteration of line search 
rl = (b - (a * xf));
ql = a * pf;
alpl = (rlt * pf) / (pft * qlt);
xl = xf + (alpl * pf);
rl1 = rl - (alpl * ql);
btl = - (rlt1 * ql) / (pft * qlt);
pl = rl1 + ( btl * pf);

% iteration of cg
qc = a * pl;
alpc = (rl1^2) / (plt *qct);
xc = xl + (alpc * pl);
rc = rl1 - (alpc * qc);
btc = (rc^2) / (rf^2);
pc = rc + (btc * pf);

 

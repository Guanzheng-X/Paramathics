
% declare all the symbols which you want to track as symbol
a = sym('A');
pi = sym('Pi');
pit = sym('PiT');
alp = sym('ALP');
ri = sym('Ri');
xi = sym('Xi');
ri1 = sym('Ri1');
bt = sym('BT');
b = sym('B');
piit = sym('PiiT');

% declare q transpose so avoid cancellation
qit = sym('QiT');
qft = sym('QfT');
qlt = sym('QlT');
qct = sym('QcT');
rlt = sym('RlT');
pft = sym('PfT');
plt = sym('PlT');
rltt = sym('RltT');

% while defining alp make sure that you use different q then qi
% because it will cancel out while we update x and r
% this code uses "qit" -> "QiT"
% iteration of cg  without fault

qi = a * pi;
alp = (ri^2) / (piit * qit);
xi = xi + (alp * pi);
ri1 = ri - (alp * qi);
bt = (ri1^2) / (ri^2);
pi = ri1 + (bt * pi);

ri = ri1;
% Ri in ri is from previous iteration


% iteration of cg when bit is flipped in a and qf gets corrupted
qf = a * pi;
alpf = (ri^2) / (pi *qft);
xf = xi + (alpf * pi);
rf = ri - (alpf * qf);
btf = (rf^2) / (ri^2);
pf = rf + (btf * pi);

% iteration of line search 
rl = (b - (a * xf));
ql = a * pf;
alpl = (rlt * pf) / (pft * qlt);
xl = xf + (alpl * pf);
rl = rf - (alpl * ql);
btl = - (rltt * ql) / (pft * qlt);
pl = rl + ( btl * pf);

% iteration of cg
qc = a * pl;
alpc = (rl^2) / (plt *qct);
xc = xl + (alpc * pl);
rc = rl - (alpc * qc);
btc = (rc^2) / (rf^2);
pc = rc + (btc * pf);

 
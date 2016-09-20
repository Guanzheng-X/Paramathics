
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
qct = sym('QcT');
qcct = sym('QccT');
pft = sym('PfT');
pct = sym('PcT');
pcct = sym('PccT');

% while defining alp make sure that you use different q then qi
% because it will cancel out while we update x and r
% this code uses "qit" -> "QiT"

% iteration of cg when bit is flipped in a and qf gets corrupted
%qf = a * pi;
alpf = (ri^2) / (pit *qft);
xf = xi + (alpf * pi);
rf = ri - (alpf * qf);
btf = (rf^2) / (ri^2);
pf = rf + (btf * pi);


% next iteratio of cg
qc = a * pf;
alpc = (rf^2) / (pft *qct);
xc = xf + (alpc * pf);
rc = rf - (alpc * qc);
btc = (rc^2) / (rf^2);
pc = rc + (btc * pf);

% final iteration of cg

qcc = a * pc;
alpcc = (rc^2) / (pct * qcct);
xcc = xc + (alpcc * pc);
rcc = rc - (alpcc * qcc);
btcc = (rcc^2) / (rc^2);
pcc = rcc + (btcc * pc);


 

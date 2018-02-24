function X=convert2au(osc)
% converts values to atomic units
h2ev = 27.21184;     %!< 1 Hartree Converts Hartree to eV
a0 = 0.529177;      %!< Bohr Radius in Angstroem

nloss=length(osc.eloss);
nt=length(osc.qtran);
emax = osc.eloss(end);
qmax = osc.qtran(end);

emax=emax/h2ev;
qmax=qmax*a0;
de=emax/nloss;
dq=qmax/nt;

ww=zeros(1,nloss);
qq=zeros(1,nt);

for i=1:nloss
    ww(i)=(i-1)*de;
end
for i=1:nt
    qq(i)=(i-1)*dq;
end

osc.A = osc.A/h2ev/h2ev;
osc.G = osc.G/h2ev;
osc.Om = osc.Om/h2ev;
osc.Ef = osc.Ef/h2ev;
osc.eloss = ww;
osc.qtran = qq;

X=osc;

end
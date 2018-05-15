function X=convert2au(osc)
% converts values to atomic units
if strcmp( osc.model,'Drude')
    osc.A = osc.A/h2ev/h2ev;
end
osc.G = osc.G/h2ev;
osc.Om = osc.Om/h2ev;
osc.Ef = osc.Ef/h2ev;
osc.eloss = osc.eloss/h2ev;
osc.qtran = osc.qtran*a0;

if isfield(osc,'u')
    osc.u = osc.u/h2ev;
end

if isfield(osc,'egap')
    osc.egap = osc.egap/h2ev;
end

X=osc;

end
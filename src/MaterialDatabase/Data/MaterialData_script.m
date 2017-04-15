clear;
% load('Elsepa.mat');
% Z2 = Z;
% for i=1:numel(Z)
%     Z2(i).n0 = Z(i).n0*10^3;
% end

MaterialData = {};

MaterialData.Au = Make_Au;

MaterialData.Mg = Make_Mg;

MaterialData.Ag = Make_Ag;

MaterialData.Si = Make_Si;

MaterialData.Al = Make_Al;

MaterialData.Nb = Make_Nb;

MaterialData.Cu = Make_Cu;

MaterialData.C = Make_C;

MaterialData.C_diamond = Make_C_diamond;

MaterialData.C_glassy = Make_C_glassy;

MaterialData.O = Make_O;

MaterialData.H = Make_H;

MaterialData.Be = Make_Be;

MaterialData.W = Make_W;

MaterialData.Ni = Make_Ni;

MaterialData.V = Make_V;

MaterialData.Ti = Make_Ti;

MaterialData.Pd = Make_Pd;

MaterialData.Ta=Make_Ta;

MaterialData.Mo=Make_Mo;

clear Au Al Si Ag i Z Z2 O Nb Cu C W Be Mg H C_diamond C_glassy Ni Ti V Pd Mo Ta
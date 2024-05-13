function [initDNS]=import_DNS(import_file)

vars_pert=readmatrix(import_file,'OutputType','char','Delimiter',{'','=',','});
values_pert=readmatrix(import_file,'OutputType','double','Delimiter',{'','=',','});

% DTMAX
Index = find(contains(vars_pert,'DTMAX'));
initDNS.(cell2mat(vars_pert(Index,1))) = values_pert(Index,2);
% Nsteps
Index = find(contains(vars_pert,'NSTEPS'));
initDNS.(cell2mat(vars_pert(Index,1))) = values_pert(Index,2);
% INTVSAVEPLANES
Index = find(contains(vars_pert,'INTVSAVEPLANES'));
initDNS.(cell2mat(vars_pert(Index,1))) = values_pert(Index,2);
% FREQ
Index = find(contains(vars_pert,'PERT_F'));
initDNS.(cell2mat(vars_pert(Index,1))) = values_pert(Index,2);
% AMP
Index = find(contains(vars_pert,'PERT_AMPL'));
initDNS.(cell2mat(vars_pert(Index,1))) = values_pert(Index,2);
initDNS.(cell2mat(vars_pert(Index,1))) = values_pert(Index,3);
% BETA
Index = find(contains(vars_pert,'PERT_BETA'));
initDNS.(cell2mat(vars_pert(Index,1))) = values_pert(Index,2);
% PERT_ZLEN
Index = find(contains(vars_pert,'PERT_ZLEN'));
initDNS.(cell2mat(vars_pert(Index,1))) = values_pert(Index,2);
% PERT_REMID
Index = find(contains(vars_pert,'PERT_REMID'));
initDNS.(cell2mat(vars_pert(Index,1))) = values_pert(Index,2);
% SPINLLEN
Index = find(contains(vars_pert,'SPINLLEN'));
initDNS.(cell2mat(vars_pert(Index,1))) = values_pert(Index,2);
% SPOUTLEN
Index = find(contains(vars_pert,'SPOUTLEN'));
initDNS.(cell2mat(vars_pert(Index,1))) = values_pert(Index,2);

end


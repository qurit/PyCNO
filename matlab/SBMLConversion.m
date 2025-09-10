%% Export model to SBML
simbiology('C:\Users\jfowler\OneDrive - UBC\Documents\Qurit\PBPK\Refactored Model\PBPK_177Lu_PSMA.sbproj')
sbr=sbioroot;
mc = sbr.Models(end);
mn = copyobj(mc);
sbmlexport(mn, 'PBPK_177Lu_PSMA.sbml')


%% need to add variants
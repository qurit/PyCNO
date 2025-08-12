%% Export model to SBML
simbiology('AlbuminModelNewImplementation.sbproj')
sbr=sbioroot;
mc = sbr.Models(end);
mn = copyobj(mc);
sbmlexport(mn, 'PSMAModel.sbml')

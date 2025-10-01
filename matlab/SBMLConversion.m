%% Export model to SBML
% simbiology('path_to_model.sbproj')
sbr=sbioroot;
mc = sbr.Models(end);
mn = copyobj(mc);
commit(mc.variants, mn); % repeat for all variants if needed by indexing varinat number i. e. commit(mc.variants(1), mn)
sbmlexport(mn, 'output.sbml')
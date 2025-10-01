%% Export model to SBML
% simbiology('AlbuminModelNewImplementation.sbproj')
sbr=sbioroot;
mc = sbr.Models(end);
mn = copyobj(mc);
commit(mc.variants, mn); % specify variant(s) number if multiple variants exist i.e. commit(mc.variants(1), mn)
sbmlexport(mn, 'output.sbml')
initCobraToolbox(false)
fileName= 'iJO1366.mat';
model= readCbModel(fileName);
model = removeRxns(model, {'BIOMASS_Ec_iJO1366_WT_53p95M'});
[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);
FBAsolution= optimizeCbModel(modelIrrev, 'max');
biomassConstraint= FBAsolution.f * 0.9;
objectiveIndex= find(modelIrrev.c== 1);
modelIrrev.lb(objectiveIndex)= biomassConstraint;
modelIrrev.c(objectiveIndex)= 0;
tempModel= modelIrrev;
outmodel = writeCbModel(modelIrrev, 'format','mat', 'fileName', 'ModelIrrevOpt90.mat');
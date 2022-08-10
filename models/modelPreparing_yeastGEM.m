initCobraToolbox(false)
fileName= 'yeast-GEM.mat';
model= readCbModel(fileName);
[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);
FBAsolution= optimizeCbModel(modelIrrev, 'max');
biomassConstraint= FBAsolution.f * 0.9;
objectiveIndex= find(modelIrrev.c== 1);
modelIrrev.lb(objectiveIndex)= biomassConstraint;
modelIrrev.c(objectiveIndex)= 0;
tempModel= modelIrrev;
outmodel = writeCbModel(modelIrrev, 'format','mat', 'fileName', 'yeast_IrrevOpt90.mat');
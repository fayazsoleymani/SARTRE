# SARTRE (ShAdow pRice - based meTabolite pRotein intEraction)
## Integration of machine learning and constraint-based modelling accurately predicts metabolite-protein interactions

SARTRE framework investigates the power of shadow prices which are calculated based on constraint-based modeling of genome-scale metabolic models(GEMs). SARTRE framework investigates the power of shadow prices which are calculated based on constraint-based modeling of genome-scale metabolic models. We can seperate the framework into five stages:

### 1. Curating the GEMs:
Two models have been used in this study. 
- [iJO1366](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3261703/): A comprehensive genome-scale reconstruction of *Escherichia coli* metabolism
- [Yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM): A consensus *S. cerevisiae* metabolic model Yeast8 and its ecosystem for comprehensively probing cellular metabolism
Curation of these mentioned GEMs have been performed by the  COBRAToolbox in MATLAB. `modelPreparing_iJO1366.m` uses `iJO1366.mat` as the input and creates `ModelIrrevOpt90.mat` and  

### 2. Calculating shadow prices

### 3. Data pre-processing

### 4. Training classifier and evaluation

### 5. Performance of SARTRE on specific tasks  





[hello](https://www.google.com/)

# SARTRE (ShAdow pRice - based meTabolite pRotein intEraction)
## Integration of machine learning and constraint-based modelling accurately predicts metabolite-protein interactions

SARTRE framework investigates the power of shadow prices which are calculated based on constraint-based modeling of genome-scale metabolic models(GEMs). SARTRE framework investigates the power of shadow prices which are calculated based on constraint-based modeling of genome-scale metabolic models. We can seperate the framework into five stages:

### 1. Curating the GEMs:
Two models have been used in this study. 
- [iJO1366](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3261703/): A comprehensive genome-scale reconstruction of *Escherichia coli* metabolism
- [Yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM): A consensus *S. cerevisiae* metabolic model Yeast8 and its ecosystem for comprehensively probing cellular metabolism
Curation of these mentioned GEMs have been performed by the  COBRAToolbox in MATLAB. `modelPreparing_iJO1366.m` uses `iJO1366.mat` as the input and creates `ModelIrrevOpt90.mat` and  `modelPreparing_yeastGEM.m` uses `Yeast-GEM.mat` and creates `yeast_IrrevOpt90.mat`.

### 2. Calculating shadow prices
In this study, we employed four gold standards integrated to GEMs as:
- [piazza](https://pubmed.ncbi.nlm.nih.gov/29307493/)
- [reznik](https://pubmed.ncbi.nlm.nih.gov/28903046/)
- [stitch_ecoli](http://stitch.embl.de/)
- [stitch_yeast](http://stitch.embl.de/)
Calculation of shadow prices can be performed in their directory by changing the working directory executing `optimization.py`. For instance, the process for piazza can be done by following commands in the terminal:
```
cd piazza
python3 optimization.py
```
The results of optimization will be saved on the `opt` directory for each reaction by its number. The process is the same for other gold standards.




### 3. Data pre-processing

### 4. Training classifier and evaluation

### 5. Performance of SARTRE on specific tasks  





[hello](https://www.google.com/)

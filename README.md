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

The optimization process has been implemented by [GUROBI](https://www.gurobi.com/) optimizer and its python interface. Make sure to be installed and have an active license for it beforehand. Calculation of shadow prices can be performed in their directory by changing the working directory followed by executing `optimization.py`. For instance, the process for `piazza` can be done by following commands in the terminal:
```
cd piazza
python3 optimization.py
```
The results of optimization will be saved on the `opt` directory for each reaction by its number. The process is the same for other gold standards. Also the results are provided in the `.zip` archive in their directory.

### 3. Data pre-processing
The framework is followed by preprocessing calculated shadow prices for downstream classifiers changing the directory to specific gold standard and executing `cleaning.py`. Also, in this stage, fingerprints of size 128 bits can be generated by executing `generate_fp.py`. make sure to have installed [RDKit](https://www.rdkit.org/) in your environment. For example for `piazza` gold standard the process can be done by executing following commands in the terminal:
```
cd piazza
python3 cleaning.py
python3 generate_fp.py
```
`generate_fp.py` prompts you to enter your email for using the [Bio.Entrez](https://biopython.org/docs/1.76/api/Bio.Entrez.html) package. The results will be saved in `met_sps_t80_replaced.csv` and `met_fp128.csv`, respectively.

### 4. Training classifier and evaluation
Preprocessing, constructing datasets, training the random forest classifier and evaluating it can be executed by running `evaluate.py` in each gold standard directory. The classifier and metrics are employed from [sklearn](https://scikit-learn.org/stable/) python library. For example, for `piazza` gold standard the process can be done by executing:
```
cd piazza
python3 evaluate.py
```
Gold standard from STITCH also get an additional argument as confidence score(150: low, 400: medium, 700: high, 900: highest) to perform this process. for example, for `stitch_ecoli` and using medium confidence score:
```
cd stitch_ecoli
python3 evaluate.py 400
```

### 5. Performance of SARTRE on specific tasks 
Additional evaluation has been provided to showcase the power of SARTRE. First, evaluation has been done by excluding shared metabolite-protein pairs, which exist in two GEMs, from the datasets, and training two separate models in remaining pairs. This can be performed by executing:
```
cd subsys_shared
python3 shared.py
```
The accuracy of two models on the test set and cosine similarity of predictions to classifiers will be displayed.
Second, we excluded metabolite-protein pairs of two subsystems of *E. coli* separately from stitch_ecoli and trained two models on remaining pairs. These subsystems are Alternate Carbon Metabolism and Cofactor and Prosthetic Group Biosynthesis. This process can be done by executing:
```
cd subsys_shared
python3 subsys.py acm
python3 subsys.py cpgb
```


***make sure to unzip `4932.protein_chemical.links.v5.0.zip archive` in both `stitch_yeast` and `subsys_shared` directories***



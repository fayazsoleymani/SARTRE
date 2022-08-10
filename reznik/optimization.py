import csv
import numpy as np
import scipy.io as sio
import gurobipy as gp
from gurobipy import GRB
import cobra
from collections import Counter
import os
import time

d = 0.2

def gpr_to_list(text):
    dic= {"or": ",",
          "(": "",
          ")": "",
          "and": ",",
          " ":""}
    for key, value in dic.items():
        text= text.replace(key, value)
    return list(set(text.split(",")))

def nparray_to_list(nparray):
    return [x[0][0] for x in nparray]

mat = sio.loadmat("ModelIrrevOpt90.mat")

S = mat['model']['S'][0][0].astype(float)
print("S: ", S.shape)
b = mat['model']['b'][0][0].reshape((1805,)).astype(float)
print("b: ", b.shape)
c = mat['model']['c'][0][0].reshape((3218,))
print("c: ", c.shape)
ub = mat['model']['ub'][0][0].reshape((3218,))
print("ub: ", ub.shape)
lb = mat['model']['lb'][0][0].reshape((3218,))
print("lb: ", lb.shape)

reactions= nparray_to_list(mat['model']['rxns'][0][0])
print("Reactions: ", len(reactions))

metabolites= nparray_to_list(mat['model']['mets'][0][0])
print("Metabolites: ", len(metabolites))

metabolite_names= nparray_to_list(mat['model']['metNames'][0][0])
genes= nparray_to_list(mat['model']['genes'][0][0])
print("Genes: ", len(genes))

gprs = [[] if len(x[0]) == 0 else gpr_to_list(x[0][0]) for x in  mat['model']['grRules'][0][0]]
print("GPRs: ", len(gprs))

article_mets= []
with open('mets.csv', 'r') as file:
    reader= csv.reader(file, delimiter= ",")
    for row in reader:
        article_mets.append(row[0])

article_genes= []
bid_to_ec_dict= dict()
with open('genes_final.csv', 'r') as file:
    reader= csv.reader(file, delimiter= ",")
    next(reader, None)
    for row in reader:
        if row[1] != '':
            article_genes.append(row[1])
            bid_to_ec_dict[row[1]]= row[0]

mets_intersection= set()
mets_mapping= []
for met in metabolites:
    if '_'.join(met.split("_")[:-1]) in article_mets:
        mets_intersection.add(met)
        mets_mapping.append([met, '_'.join(met.split("_")[:-1])])


with open('mets_final.csv', "w", newline='', encoding= "UTF8") as file:
    writer= csv.writer(file, delimiter = ',')
    writer.writerow(['ModelMet', 'ArticleMet'])
    writer.writerows(mets_mapping)

target_reaction_indeces = []
for i in range(len(reactions)):
    if len(gprs[i]) != 0:
        target_reaction_indeces.append(i)


target_metabolites_indeces = []
for i in range(len(metabolites)):
    if metabolites[i] in mets_intersection:
        target_metabolites_indeces.append(i)


target_metabolites_names= []
for index in target_metabolites_indeces:
    target_metabolites_names.append([metabolites[index]])
len(target_metabolites_names)
with open('target_metabolites.tsv', "w", newline='', encoding= "UTF8") as file:
    writer= csv.writer(file, delimiter = '\t')
    writer.writerows(target_metabolites_names)

target_reaction_names= []
for index in target_reaction_indeces:
    target_reaction_names.append([reactions[index]])
len(target_reaction_names)
with open('target_reactions.tsv', "w", newline='', encoding= "UTF8") as file:
    writer= csv.writer(file, delimiter = '\t')
    writer.writerows(target_reaction_names)

model_genes= []
for gpr in gprs:
    model_genes.extend(gpr)
model_genes= set(model_genes)
len(model_genes)

target_genes= list(model_genes.intersection(set(article_genes)))
print("Intersection of genes between model and gold standard:", len(target_genes))

target_genes_names= []
for g in target_genes:
    target_genes_names.append([g])
with open('target_genes.tsv', "w", newline='', encoding= "UTF8") as file:
    writer= csv.writer(file, delimiter = '\t')
    writer.writerows(target_genes_names)

gene_in_rxns= []
for gene in target_genes:
    temp= [bid_to_ec_dict[gene]]
    for index, rxn in enumerate(reactions):
        if gene in gprs[index]:
            temp.append(1)
        else:
            temp.append(0)
    gene_in_rxns.append(temp)
with open('gene_in_rxns.csv', "w", newline='', encoding= "UTF8") as file:
    writer= csv.writer(file, delimiter = ',')
    writer.writerows(gene_in_rxns)

print("The optimization will be performed to calculate shadow prices in {} reactions for {} metabolites".format(len(target_reaction_indeces), len(target_metabolites_indeces)))


os.mkdir('./opt')

for objective_index in target_reaction_indeces:
    start_time= time.time()
    c[objective_index] = 1

    with gp.Env(empty=True) as env:
        env.setParam('OutputFlag', 0)
        env.start()
        with gp.Model(env=env) as model:
            v= model.addMVar(c.shape, vtype= GRB.CONTINUOUS, lb= lb, ub= ub)
            steady_state= model.addMConstr(S, v, '=', b)
            model.setObjective(c @ v, GRB.MAXIMIZE)
            model.optimize()
            objective_value= model.ObjVal
            shadow_prices= steady_state.getAttr(GRB.Attr.Pi)
            g_plus= steady_state.getAttr(GRB.Attr.SARHSUp)
            g_minus= steady_state.getAttr(GRB.Attr.SARHSLow)


    decremental_shadow_prices = dict()
    for i in target_metabolites_indeces:
        if g_minus[i] == 0:
            decremental_shadow_price = "Inf"
        else:
            perturbation = d * g_minus[i]
            b[i] = perturbation
            with gp.Env(empty=True) as env:
                env.setParam('OutputFlag', 0)
                env.start()
                with gp.Model(env=env) as model:
                    v= model.addMVar(c.shape , vtype= GRB.CONTINUOUS, lb= lb, ub= ub)
                    steady_state= model.addMConstr(S, v, '=', b)
                    model.setObjective(c @ v, GRB.MAXIMIZE)
                    model.optimize()
                    if model.status == GRB.OPTIMAL:
                        decremental_shadow_price= (model.ObjVal - objective_value) / perturbation
            b[i] = 0
        decremental_shadow_prices[i] = decremental_shadow_price


    incremental_shadow_prices = dict()
    for i in target_metabolites_indeces:
        if g_plus[i] == 0:
            incremental_shadow_price = "Inf"
        else:
            perturbation = d * g_plus[i]
            b[i] = perturbation
            with gp.Env(empty=True) as env:
                env.setParam('OutputFlag', 0)
                env.start()
                with gp.Model(env=env) as model:
                    v= model.addMVar(c.shape , vtype= GRB.CONTINUOUS, lb= lb, ub= ub)
                    steady_state= model.addMConstr(S, v, '=', b)
                    model.setObjective(c @ v, GRB.MAXIMIZE)
                    model.optimize()
                    if model.status == GRB.OPTIMAL:
                        incremental_shadow_price= (model.ObjVal - objective_value) / perturbation
            b[i] = 0

        incremental_shadow_prices[i] = incremental_shadow_price

    header= ['Reaction', 'Metabolite', 'ShadowPrice', 'G-', 'G+', 'DecrementalShadowPrice', 'IncrementalShadowPrice']
    rows= []
    for i in target_metabolites_indeces:
        rows.append([reactions[objective_index], metabolites[i], shadow_prices[i],
                     g_minus[i], g_plus[i], decremental_shadow_prices[i], incremental_shadow_prices[i]])
    with open(os.path.join("opt", str(objective_index) + '.tsv'),
              "w", newline='', encoding= "UTF8") as file:
        writer= csv.writer(file, delimiter = '\t')
        writer.writerow(header)
        writer.writerows(rows)


    c[objective_index] = 0


    end_time= time.time()
    print(objective_index, ':', end_time- start_time)

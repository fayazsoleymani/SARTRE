import csv
import numpy as np
import scipy.io as sio
import gurobipy as gp
from gurobipy import GRB
import cobra
from collections import Counter
import os

d = 0.2

def gpr_to_list(text):
    dic= {"or": ",",
          "(": "",
          ")": "",
          "and": ",",
          " ":""}
    for key, value in dic.items():
        text= text.replace(key, value)
    return text.split(",")

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

uniprot_to_bid_dict= dict()
b_id_to_uniprot_dict=dict()
uniprot_ids= set()
b_ids= set()
with open("mapped_genes_mmc5.csv", "r") as file:
    reader= csv.reader(file, delimiter= ",")
    next(reader, None)
    for row in reader:
        uniprot_to_bid_dict[row[0]] = row[1]
        b_id_to_uniprot_dict[row[1]]= row[0]
        uniprot_ids.add(row[0])
        b_ids.add(row[1])


article_to_model_mets_dict= dict()
model_mets_to_article_dict= dict()
article_mets= set()
model_mets= set()
with open("mapped_metabolites.tsv", "r") as file:
    reader= csv.reader(file, delimiter= "\t")
    next(reader, None)
    for row in reader:
        article_to_model_mets_dict[row[0]] = row[1]
        model_mets_to_article_dict[row[1]]= row[0]
        article_mets.add(row[0])
        model_mets.add(row[1])


target_reaction_indeces = []
for i in range(len(reactions)):
    if len(gprs[i]) != 0:
        target_reaction_indeces.append(i)


target_reaction_names= []
for index in target_reaction_indeces:
    target_reaction_names.append([reactions[index]])
len(target_reaction_names)
with open('target_reactions.tsv', "w", newline='', encoding= "UTF8") as file:
    writer= csv.writer(file, delimiter = '\t')
    writer.writerows(target_reaction_names)

target_metabolites_indeces = []
mets_final= []
for i in range(len(metabolites)):
    if metabolite_names[i] in model_mets:
        target_metabolites_indeces.append(i)
        mets_final.append([metabolites[i], model_mets_to_article_dict[metabolite_names[i]]])


target_metabolites_names= []
for index in target_metabolites_indeces:
    target_metabolites_names.append([metabolites[index]])
len(target_metabolites_names)
with open('target_metabolites.tsv', "w", newline='', encoding= "UTF8") as file:
    writer= csv.writer(file, delimiter = '\t')
    writer.writerows(target_metabolites_names)

with open('mets_final.csv', 'w', newline='', encoding= "UTF8") as file:
    writer= csv.writer(file, delimiter = ',')
    writer.writerow(['ModelMet', 'ArticleMet'])
    writer.writerows(mets_final)

target_genes= list(set(genes).intersection(set(b_ids)))
print("Intersection of genes between model and gold standard:", len(target_genes))


gene_in_rxns= []
for gene in target_genes:
    temp= [b_id_to_uniprot_dict[gene]]
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
        rows.append([reactions[objective_index], metabolites[i], shadow_prices[i], g_minus[i], g_plus[i], decremental_shadow_prices[i], incremental_shadow_prices[i]])
    with open(os.path.join("opt", str(objective_index) + '.tsv'), "w", newline='', encoding= "UTF8") as file:
        writer= csv.writer(file, delimiter = '\t')
        writer.writerow(header)
        writer.writerows(rows)

    c[objective_index] = 0

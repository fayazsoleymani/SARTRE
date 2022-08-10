import csv
import numpy as np
import scipy.io as sio
import gurobipy as gp
from gurobipy import GRB
from collections import Counter
import os
import time

d = 0.2

def gpr_to_list(text):
    dic= {"or": "|",
          "(": "",
          ")": "",
          "and": "|",
          " ":""}
    for key, value in dic.items():
        text= text.replace(key, value)
    return list(set(text.split("|")))

def nparray_to_list(nparray):
    return [x[0][0] for x in nparray]

mat = sio.loadmat('./yeast_IrrevOpt90.mat')

S = mat['model']['S'][0][0].astype(float)
print("S: ", S.shape)
b = mat['model']['b'][0][0].reshape((-1)).astype(float)
print("b: ", b.shape)
c = mat['model']['c'][0][0].reshape((-1))
print("c: ", c.shape)
ub = mat['model']['ub'][0][0].reshape((-1))
print("ub: ", ub.shape)
lb = mat['model']['lb'][0][0].reshape((-1))
print("lb: ", lb.shape)

reactions= nparray_to_list(mat['model']['rxns'][0][0])
print("Reactions: ", len(reactions))

metabolites= nparray_to_list(mat['model']['mets'][0][0])
print("Metabolites: ", len(metabolites))

metabolite_names= nparray_to_list(mat['model']['metNames'][0][0])

genes= nparray_to_list(mat['model']['genes'][0][0])
print("Genes: ", len(genes))

gprs= []
with open("gpr_rules.tsv", "r") as file:
    reader= csv.reader(file, delimiter= "\t")
    for row in reader:
        if row[1] == '':
            gprs.append([])
        else:
            gprs.append(gpr_to_list(row[1]))
print("GPRs: ", len(gprs))

stitch_mets= set()
stitch_genes= set()
with open("4932.protein_chemical.links.v5.0.tsv", "r") as file:
    reader= csv.reader(file, delimiter= "\t")
    next(reader, None)
    for line in reader:
        met= line[0]
        if "CIDm" in met:
            met= met.replace("CIDm", "")
        elif "CIDs" in met:
            met= met.replace("CIDs", "")
        met= str(int(met))
        stitch_mets.add(met)
        gene= line[1].split(".")[1]
        stitch_genes.add(gene)

model_cids= set()
mets_to_cid= dict()
with open('mets_final.csv', 'r') as file:
    reader= csv.reader(file, delimiter= ",")
    next(reader, None)
    for row in reader:
        if row[2] != '':
            mets_to_cid[row[0]]= row[2]
            model_cids.add(row[2])

mets_intersection= model_cids.intersection(stitch_mets)

target_reaction_indeces = []
for i in range(len(reactions)):
    if len(gprs[i]) != 0:
        target_reaction_indeces.append(i)

target_metabolites_indeces = []
for i in range(len(metabolites)):
    if metabolites[i] in mets_to_cid:
        if mets_to_cid[metabolites[i]] in mets_intersection:
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

target_genes= list(set(genes).intersection(stitch_genes))
print("Intersection of genes between model and gold standard:", len(target_genes))

target_genes_names= []
for g in target_genes:
    target_genes_names.append([g])
with open('target_genes.tsv', "w", newline='', encoding= "UTF8") as file:
    writer= csv.writer(file, delimiter = '\t')
    writer.writerows(target_genes_names)

gene_name_to_id= dict()
with open("gene_names_to_id.csv", "r") as file:
    reader= csv.reader(file, delimiter= ",")
    next(reader, None)
    for row in reader:
        gene_name_to_id[row[0]]= row[1]
len(gene_name_to_id)

gprs_based_on_id= []
for entry in gprs:
    temp= []
    for gene in entry:
        temp.append(gene_name_to_id[gene])
    gprs_based_on_id.append(temp)
len(gprs_based_on_id)

gene_in_rxns= []
for gene in target_genes:
    temp= [gene]
    for index, rxn in enumerate(reactions):
        if gene in gprs_based_on_id[index]:
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

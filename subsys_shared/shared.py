import csv
import numpy as np
from collections import Counter
from imblearn.under_sampling import RandomUnderSampler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics.pairwise import cosine_similarity


def make_float(a):
    temp= []
    for i in a:
        temp_row= []
        for j in i:
            if j == '':
                temp_row.append(np.nan)
            else:
                temp_row.append(float(j))
        temp.append(temp_row)
    return np.array(temp)


def read_file(file_name):
    data= []
    with open(file_name, 'r') as file:
        reader= csv.reader(file, delimiter= ',')
        for row in reader:
            data.append(row)
    data= np.array(data)
    names= data[:, 0].reshape(-1,1)
    features= data[:, 1:].astype('str')

    return names, make_float(features)


def make_unique_without_zeros(features, decimal_points):
    unique_features_set= []
    for j in range(features.shape[1]):
        feature_set= []
        for i in range(features.shape[0]):
            feature_set.append(round(features[i][j], decimal_points))
        if all(x == 0 for x in feature_set):
            continue
        if feature_set not in unique_features_set:
            unique_features_set.append(feature_set)
    return np.array(unique_features_set).transpose()


def read_stitch_interaction_file(file_name, met_names, gene_names, confidence_score):
    interactions= dict()
    with open(file_name, 'r') as file:
            reader= csv.reader(file, delimiter= "\t")
            next(reader, None)
            for row in reader:
                met= row[0]
                if "CIDm" in met:
                    met= met.replace("CIDm", "")
                elif "CIDs" in met:
                    met= met.replace("CIDs", "")
                met= str(int(met))
                gene= row[1].split(".")[1]
                if met in met_names and gene in gene_names:
                    if int(row[2]) >= confidence_score:
                        interaction= 1
                    else:
                        interaction= 0
                    interactions[(met, gene)]= interaction
    return interactions


def delete_redundant_entries(names, features):
    unique_names= []
    for name in names:
        if name[0] not in unique_names:
            unique_names.append(name[0])
    new_names=[]
    new_features= []
    for name in unique_names:
        temp= []
        for index, entry in enumerate(names.reshape(-1)):
            if name == entry:
                temp.append(features[index])
        summation= temp[0]
        count= len(temp)
        for x in temp[1:]:
            summation+=x
        temp= summation/count

        new_names.append(name)
        new_features.append(temp)
    return np.array(new_names).reshape(-1, 1), np.array(new_features)



def read_specified_file(filename, target, decimal_point):
    names, features= read_file(filename)
    features= make_unique_without_zeros(features, decimal_point)
    names, features= delete_redundant_entries(names, features)

    new_names= []
    new_features= []
    for entry in target:
        for index, name in enumerate(names):
            if entry == name[0]:
                new_names.append(entry)
                new_features.append(features[index])
    return np.array(new_names).reshape(-1, 1), np.array(new_features)


def make_dataset(met_names, met_features, gene_names, gene_features, gold_standard, confidence_score):

    met_names= met_names.transpose()[0].tolist()
    gene_names= gene_names.transpose()[0].tolist()

    interactions= dict()

    if gold_standard == 'stitch_ecoli':
        interactions= read_stitch_interaction_file('./511145.protein_chemical.links.v5.0.tsv',
                                                    met_names, gene_names, confidence_score)

    elif gold_standard == 'stitch_yeast':
        interactions= read_stitch_interaction_file('./4932.protein_chemical.links.v5.0.tsv',
                                                   met_names, gene_names, confidence_score)



    dataset= np.zeros(( met_features.shape[1]+ gene_features.shape[1]+ 1))
    for i, m in enumerate(met_names):
        temp = []
        for j, p in enumerate(gene_names):
            if (m, p) in interactions:
                interaction= interactions[(m, p)]
            else:
                interaction= 0
            l= met_features[i].tolist()
            l.extend(gene_features[j])
            l.append(interaction)

            temp.append(l)
        dataset= np.vstack((dataset, np.array(temp)))

    dataset= np.delete(dataset, 0, axis= 0)
    print("Dataset shape:\t", dataset.shape)

    X= dataset[:, :-1]
    y= dataset[:, -1]

    print("X shape:\t", X.shape, "\ny shape:\t", y.shape)
    return X, y

def make_dataset_exception(met_names, met_features, gene_names, gene_features, gold_standard, confidence_score, exceptions):

    met_names= met_names.transpose()[0].tolist()
    gene_names= gene_names.transpose()[0].tolist()

    interactions= dict()


    if gold_standard == 'stitch_ecoli':
        interactions= read_stitch_interaction_file('./511145.protein_chemical.links.v5.0.tsv',
                                                    met_names, gene_names, confidence_score)

    elif gold_standard == 'stitch_yeast':
        interactions= read_stitch_interaction_file('./4932.protein_chemical.links.v5.0.tsv',
                                                   met_names, gene_names, confidence_score)

    else:
        print("gold standard name is not valid")
        return

    print("Number of entries in interaction file:\t",  len(interactions))


    dataset= np.zeros(( met_features.shape[1]+ gene_features.shape[1]+ 1))
    for i, m in enumerate(met_names):
        temp = []
        for j, p in enumerate(gene_names):
            if (m, p) not in exceptions:
                if (m, p) in interactions:
                    interaction= interactions[(m, p)]
                else:
                    interaction= 0
                l= met_features[i].tolist()
                l.extend(gene_features[j])
                l.append(interaction)

                temp.append(l)
        dataset= np.vstack((dataset, np.array(temp)))

    dataset= np.delete(dataset, 0, axis= 0)
    print("Dataset shape:\t", dataset.shape)

    X= dataset[:, :-1]
    y= dataset[:, -1]

    print("X shape:\t", X.shape, "\ny shape:\t", y.shape)
    return X, y


def make_dataset_exact(met_names, met_features, gene_names, gene_features, gold_standard, confidence_score, exacts):

    met_names= met_names.transpose()[0].tolist()
    gene_names= gene_names.transpose()[0].tolist()

    interactions= dict()



    if gold_standard == 'stitch_ecoli':
        interactions= read_stitch_interaction_file('./511145.protein_chemical.links.v5.0.tsv',
                                                    met_names, gene_names, confidence_score)

    elif gold_standard == 'stitch_yeast':
        interactions= read_stitch_interaction_file('./4932.protein_chemical.links.v5.0.tsv',
                                                   met_names, gene_names, confidence_score)

    else:
        print("gold standard name is not valid")
        return

    print("Number of entries in interaction file:\t",  len(interactions))


    dataset= np.zeros(( met_features.shape[1]+ gene_features.shape[1]+ 1))
    for i, m in enumerate(met_names):
        temp = []
        flag= False
        for j, p in enumerate(gene_names):
            if (m, p) in exacts:
                flag= True
                if (m, p) in interactions:
                    interaction= interactions[(m, p)]
                else:
                    interaction= 0
                l= met_features[i].tolist()
                l.extend(gene_features[j])
                l.append(interaction)

                temp.append(l)
        if flag:
            dataset= np.vstack((dataset, np.array(temp)))

    dataset= np.delete(dataset, 0, axis= 0)
    print("Dataset shape:\t", dataset.shape)

    X= dataset[:, :-1]
    y= dataset[:, -1]

    print("X shape:\t", X.shape, "\ny shape:\t", y.shape)
    return X, y

stitch_cid_intersection= ['1117', '312', '5793', '888']

ec_intersection= []
ecoli_gene_intersection= []
yeast_gene_intersection= []
with open('./stitch_ec_intersections.csv', 'r') as file:
    reader= csv.reader(file, delimiter= ',')
    next(reader, None)
    for row in reader:
        ec_intersection.append(row[0])
        ecoli_gene_intersection.append(row[1])
        yeast_gene_intersection.append(row[2])

ecoli_interactions= read_stitch_interaction_file('./511145.protein_chemical.links.v5.0.tsv', stitch_cid_intersection, ecoli_gene_intersection, 400)
yeast_interactions= read_stitch_interaction_file('./4932.protein_chemical.links.v5.0.tsv', stitch_cid_intersection, yeast_gene_intersection, 400)

valid_ecoli_pairs= []
valid_yeast_pairs= []
for cid in stitch_cid_intersection:
    for ec, ecoli_gene, yeast_gene in zip(ec_intersection, ecoli_gene_intersection, yeast_gene_intersection):
        ecoli_interaction= ecoli_interactions[(cid, ecoli_gene)] if (cid, ecoli_gene) in ecoli_interactions else 0
        yeast_interaction= yeast_interactions[(cid, yeast_gene)] if (cid, yeast_gene) in yeast_interactions else 0
        if ecoli_interaction == yeast_interaction:
            valid_ecoli_pairs.append((cid, ecoli_gene))
            valid_yeast_pairs.append((cid, yeast_gene))
len(valid_ecoli_pairs)

ecoli_exceptions= [(m, p) for m in stitch_cid_intersection for p in ecoli_gene_intersection if (m,p) in valid_ecoli_pairs]


ecoli_cid_filename= './ecoli/met_sps_t80_replaced.csv'
ecoli_gene_filename= './ecoli/gene_in_rxns.csv'
ecoli_all_cid_names, ecoli_all_cid_features = read_file(ecoli_cid_filename)
ecoli_all_cid_features= make_unique_without_zeros(ecoli_all_cid_features, 2)
ecoli_all_cid_names, ecoli_all_cid_features= delete_redundant_entries(ecoli_all_cid_names, ecoli_all_cid_features)
ecoli_all_gene_names, ecoli_all_gene_features =read_file(ecoli_gene_filename)
ecoli_all_gene_features= make_unique_without_zeros(ecoli_all_gene_features, 0)
ecoli_intersection_cid_names, ecoli_intersection_cid_features= read_specified_file(ecoli_cid_filename, stitch_cid_intersection, 2)
ecoli_intersection_gene_names, ecoli_intersection_gene_features= read_specified_file(ecoli_gene_filename, ecoli_gene_intersection, 0)

print("Ecoli Test:")
ecoli_X_test, ecoli_y_test= make_dataset_exact(ecoli_intersection_cid_names, ecoli_intersection_cid_features,
                                               ecoli_intersection_gene_names, ecoli_intersection_gene_features,
                                               'stitch_ecoli', 400, ecoli_exceptions)

print("Ecoli_train:")
ecoli_X_train, ecoli_y_train= make_dataset_exception(ecoli_all_cid_names, ecoli_all_cid_features,
                                                     ecoli_all_gene_names, ecoli_all_gene_features,
                                                     'stitch_ecoli', 400, ecoli_exceptions)


rus= RandomUnderSampler(random_state=42)
ecoli_X_res, ecoli_y_res= rus.fit_resample(ecoli_X_train, ecoli_y_train)
print('Original Set:', Counter(ecoli_y_train))
print('Resampled Set:', Counter(ecoli_y_res))
clf= RandomForestClassifier(n_estimators= 100)
clf.fit(ecoli_X_res, ecoli_y_res)
ecoli_y_pred= clf.predict(ecoli_X_test)
print("Ecoli accuracy on shared test set:", accuracy_score(ecoli_y_test, ecoli_y_pred))


yeast_exceptions= [(m, p) for m in stitch_cid_intersection for p in yeast_gene_intersection if (m,p) in valid_yeast_pairs]

yeast_cid_filename= './yeast/met_sps_t80_replaced.csv'
yeast_gene_filename= './yeast/gene_in_rxns.csv'



yeast_all_cid_names, yeast_all_cid_features = read_file(yeast_cid_filename)
yeast_all_cid_features= make_unique_without_zeros(yeast_all_cid_features, 2)
yeast_all_cid_names, yeast_all_cid_features= delete_redundant_entries(yeast_all_cid_names, yeast_all_cid_features)

yeast_all_gene_names, yeast_all_gene_features =read_file(yeast_gene_filename)
yeast_all_gene_features= make_unique_without_zeros(yeast_all_gene_features, 0)

yeast_intersection_cid_names, yeast_intersection_cid_features= read_specified_file(yeast_cid_filename, stitch_cid_intersection, 2)
yeast_intersection_gene_names, yeast_intersection_gene_features= read_specified_file(yeast_gene_filename, yeast_gene_intersection, 0)

print("Yeast test set:")
yeast_X_test, yeast_y_test= make_dataset_exact(yeast_intersection_cid_names, yeast_intersection_cid_features,
                                                     yeast_intersection_gene_names, yeast_intersection_gene_features,
                                                     'stitch_yeast', 400, yeast_exceptions)
print("Yeast train set:")
yeast_X_train, yeast_y_train= make_dataset_exception(yeast_all_cid_names, yeast_all_cid_features,
                                                     yeast_all_gene_names, yeast_all_gene_features,
                                                     'stitch_yeast', 400, yeast_exceptions)

rus= RandomUnderSampler(random_state= 42)
yeast_X_res, yeast_y_res= rus.fit_resample(yeast_X_train, yeast_y_train)
print('Original Set:', Counter(yeast_y_train))
print('Resampled Set:', Counter(yeast_y_res))
yeast_clf= RandomForestClassifier(n_estimators= 100)
yeast_clf.fit(yeast_X_res, yeast_y_res)
yeast_y_pred= yeast_clf.predict(yeast_X_test)
print("Yeast accuracy on shared test set:", accuracy_score(yeast_y_test, yeast_y_pred))

print("Cosine between predictions of shared test set of Ecoli and Yeast:", cosine_similarity(np.array([ecoli_y_pred]), np.array([yeast_y_pred]))[0][0])

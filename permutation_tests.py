import csv
import numpy as np
from collections import Counter
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_validate
from imblearn.under_sampling import RandomUnderSampler
from sklearn.ensemble import RandomForestClassifier
import sys
from sklearn.metrics import roc_auc_score, f1_score, accuracy_score

permutation= sys.argv[1]

metabolites_features_filename= './met_sps_t80_replaced.csv'
protein_features_filename= './gene_in_rxns.csv'


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


def make_dataset(met_names, met_features, gene_names, gene_features):

    met_names= met_names.transpose()[0].tolist()
    gene_names= gene_names.transpose()[0].tolist()

    interactions= dict()


    with open('./interactions_mmc5.csv', 'r') as file:
        reader= csv.reader(file, delimiter= ',')
        next(reader, None)
        for row in reader:
            interactions[(row[1], row[0])] = int(row[2])

    print("Number of entries in interaction file:\t",  len(interactions))


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

def feature_permutate(X, y):
    positive_indeces= []
    negative_indeces= []
    for index, label in enumerate(y):
        if label == 0:
            negative_indeces.append(index)
        elif label == 1:
            positive_indeces.append(index)
    #print(negative_indeces)
    #print(positive_indeces)
    X_negative, y_negative= X[negative_indeces], y[negative_indeces]
    X_positive, y_positive= X[positive_indeces], y[positive_indeces]
    np.random.seed= 0
    np.transpose(np.random.shuffle(np.transpose(X_negative)))
    np.random.seed= 1
    np.transpose(np.random.shuffle(np.transpose(X_positive)))

    X= np.vstack((X_negative, X_positive))
    y= np.hstack((y_negative, y_positive))

    return X, y

metabolites_names, metabolites_features= read_file(metabolites_features_filename)
print("Number of metabolites:", metabolites_features.shape[0], '\t',"Lenght of features:" , metabolites_features.shape[1])

metabolites_features= make_unique_without_zeros(metabolites_features, 2)
metabolites_names, metabolites_features= delete_redundant_entries(metabolites_names, metabolites_features)
print("Number of metabolites after preprocessing:",metabolites_features.shape[0], '\t',"Lenght of features after preprocessing:" , metabolites_features.shape[1])

protein_names, protein_features= read_file(protein_features_filename)
print("Number of proteins:", protein_features.shape[0], '\t',"Lenght of features:", protein_features.shape[1])

protein_features= make_unique_without_zeros(protein_features, 0)
print("Number of proteins after preprocessing:", protein_features.shape[0], '\t',"Lenght of features after preprocessing:", protein_features.shape[1])

X, y= make_dataset(metabolites_names, metabolites_features,
                   protein_names, protein_features)



org_aucs= []
y_true= np.array([])
y_pred= np.array([])

kf= KFold(n_splits= 10, shuffle= True)

for train_index, test_index in kf.split(X):
    X_train , X_test = X[train_index], X[test_index]
    y_train , y_test = y[train_index] , y[test_index]

    clf= RandomForestClassifier(n_estimators= 100)
    rus= RandomUnderSampler(random_state=42)

    X_train_res, y_train_res= rus.fit_resample(X_train, y_train)

    print('Original Set:', Counter(y_train))
    print('Resampled Set:', Counter(y_train_res))

    clf.fit(X_train_res, y_train_res)

    org_aucs.append(roc_auc_score(y_test, clf.predict(X_test), average='macro'))

    y_pred= np.hstack((y_pred, clf.predict(X_test)))
    y_true= np.hstack((y_true, y_test))

print("Original AUC:", roc_auc_score(y_true, y_pred, average='macro'))

permutation_aucs= []
y_true= np.array([])
y_pred= np.array([])

kf= KFold(n_splits= 10, shuffle= True)

for train_index, test_index in kf.split(X):
    X_train , X_test = X[train_index], X[test_index]
    y_train , y_test = y[train_index] , y[test_index]

    clf= RandomForestClassifier(n_estimators= 100)
    rus= RandomUnderSampler(random_state=42)

    if permutation == 'feature_permutation':
        X_train, y_train= feature_permutate(X_train, y_train)       #feature permutation

    X_train_res, y_train_res= rus.fit_resample(X_train, y_train)

    if permutation == 'label_permutation':
        np.random.shuffle(y_train_res)      #label permutation


    print('Original Set:', Counter(y_train))
    print('Resampled Set:', Counter(y_train_res))

    clf.fit(X_train_res, y_train_res)

    permutation_aucs.append(roc_auc_score(y_test, clf.predict(X_test), average='macro'))

    y_pred= np.hstack((y_pred, clf.predict(X_test)))
    y_true= np.hstack((y_true, y_test))

print("Permutated AUC:", roc_auc_score(y_true, y_pred, average='macro'))

print("P-value:\t", stats.ttest_ind(org_aucs, permutation_aucs, alternative='greater'))

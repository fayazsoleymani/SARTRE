import csv
import os


root='./opt'
target_metabolites_file= './target_metabolites.tsv'
target_reactions_file= './target_reactions.tsv'
mets_mapping_file= './mets_final.csv'
mapped_index= 1
final_output_file= './met_sps_t80_replaced.csv'



def delete_inf(entry):
    if entry[3] != "Inf" and entry[4] != "Inf":
        cleaned= [entry[0], entry[1], float(entry[2]), float(entry[3]), float(entry[4])]
    else:
        cleaned= [entry[0], entry[1], None]
    return cleaned


tolerance = 1e-05
decimal_points = 5
def delete_two_sided(entry):
    if entry[2] == None:
        return entry
    elif abs(entry[3] - entry[4]) < tolerance and round(entry[2], decimal_points) == round(entry[3], decimal_points) and round(entry[2], decimal_points) == round(entry[4], decimal_points):
        return [entry[0], entry[1], round(entry[2], decimal_points)]
    else:
        return [entry[0], entry[1], None]


def read_file(root, index):
    temp_list = []
    with open(os.path.join(root, str(index) + '.tsv'), 'r') as file:
        reader = csv.reader(file, delimiter = '\t')
        next(reader, None)
        for row in reader:
            temp_list.append([row[0], row[1], row[2], row[5], row[6]])
    cleaned_list= []
    for entry in temp_list:
        cleaned_list.append(delete_two_sided(delete_inf(entry)))
    return cleaned_list

def read_file_without_validating(root, index):
    temp_list = []
    with open(os.path.join(root, str(index) + '.tsv'), 'r') as file:
        reader = csv.reader(file, delimiter = '\t')
        next(reader, None)
        for row in reader:
            temp_list.append([row[0], row[1], round(float(row[2]), 5)])
    return temp_list

files = os.listdir(root)
print("Number of files in the root", len(files))
reaction_indeces = set()
for file in files:
    reaction_indeces.add(int(file.split(".")[0]))
reaction_indeces= list(sorted(reaction_indeces))

opt_data= []
for index in reaction_indeces:
    opt_data.extend(read_file(root, index))

target_metabolites = []
with open(target_metabolites_file, 'r') as file:
    reader = csv.reader(file, delimiter = '\t')
    for row in reader:
        target_metabolites.append(row[0])

target_reactions = []
with open(target_reactions_file, 'r') as file:
    reader = csv.reader(file, delimiter = '\t')
    for row in reader:
        target_reactions.append(row[0])
len(target_reactions)

opt_data_dict= dict()
for entry in opt_data:
    opt_data_dict[(entry[0], entry[1])]= entry[2]

met_sps_all= dict()
for entry in target_metabolites:
    met_sps_all[entry]= []
for met in target_metabolites:
    for rxn in target_reactions:
        met_sps_all[met].append(opt_data_dict[(rxn, met)])

met_sps_nonempty_count= dict()
for key, value in met_sps_all.items():
    count= 0
    for sp in value:
        if sp != None:
            count += 1
    met_sps_nonempty_count[key] = count

threshold= round(0.80*len(target_reactions))

valuable_mets= []
for key, value in met_sps_nonempty_count.items():
    if value> threshold:
        valuable_mets.append(key)

opt_data_unvalidated= []
for index in reaction_indeces:
    opt_data_unvalidated.extend(read_file_without_validating(root, index))

opt_data_unvalidated_dict= dict()
for entry in opt_data_unvalidated:
    opt_data_unvalidated_dict[(entry[0], entry[1])]= entry[2]

met_sps= dict()
for entry in valuable_mets:
    met_sps[entry]= []
for met in valuable_mets:
    for rxn in target_reactions:
        if opt_data_dict[(rxn, met)] != None:
            met_sps[met].append(opt_data_dict[(rxn, met)])
        else:
            met_sps[met].append(opt_data_unvalidated_dict[(rxn, met)])

mets_mapping= dict()
with open(mets_mapping_file, 'r') as file:
    reader= csv.reader(file, delimiter= ",")
    next(reader, None)
    for row in reader:
        if row[mapped_index] != '':
            mets_mapping[row[0]]= row[mapped_index]

rows= []
for key, value in met_sps.items():
    temp= [mets_mapping[key]]
    temp.extend(value)
    rows.append(temp)

with open(final_output_file, 'w', newline='', encoding= "UTF8") as file:
    writer= csv.writer(file, delimiter = ',')
    writer.writerows(rows)

print("Number of final validated metabolites:", len(rows))

import csv
from rdkit import Chem
from Bio import Entrez

Entrez.email= input("Enter your Entrez Email:")

base_mets_filename= './met_sps_t80_replaced.csv'
fingerprint_size= 128
cid_replacement= True
cid_replacement_filename= './mets_cid.csv'
cid_col= 1
output_file= './met_fp' + str(fingerprint_size) + '.csv'

def fetch_smiles_from_cid(cid):
    record = Entrez.read(Entrez.esummary(db="pccompound", id=cid))
    smiles= record[0]['CanonicalSmiles']
    return smiles


def retrieve_fp_from_smiles(smiles, fp_len):
    fp= Chem.RDKFingerprint(Chem.MolFromSmiles(smiles), fpSize= fp_len).ToBitString()
    return fp


def split_to_char(word):
    return [char for char in word]

base_mets= []
with open(base_mets_filename, 'r') as file:
    reader= csv.reader(file, delimiter= ',')
    for row in reader:
        base_mets.append(row[0])
print("Number of metabolites:",len(base_mets))

if cid_replacement:
    mets_cid_dict= dict()
    with open(cid_replacement_filename, 'r') as file:
        reader= csv.reader(file, delimiter= ',')
        next(reader, None)
        for row in reader:
            mets_cid_dict[row[0]]= row[cid_col]

mets_fp= []
for met in base_mets:
    cid= mets_cid_dict[met] if cid_replacement else met
    smiles= fetch_smiles_from_cid(cid)
    fp= retrieve_fp_from_smiles(smiles, fingerprint_size)
    mets_fp.append([met] + split_to_char(fp))


with open(output_file, 'w', newline='', encoding= "UTF8") as file:
    writer= csv.writer(file, delimiter = ',')
    writer.writerows(mets_fp)

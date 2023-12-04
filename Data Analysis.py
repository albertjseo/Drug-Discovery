import pandas as pd
import numpy as np
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

# PART 2: Exploratory Data Analyis
newDF = pd.read_csv('Bioactivity_preprocessed_data.csv')

#collect Lipinski descriptors which highlights the druglikeness and pharmacokinetic profiles. 
#inspired by: https://codeocean.com/explore/capsules?query=tag:data-curation
def lipinski(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem) 
        moldata.append(mol)
       
    baseData = np.arange(1,1)
    i=0  
    for mol in moldata:        
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
           
        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])   
    
        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i = i+1      
    
    columnNames = ["MW","LogP","NumHDonors","NumHAcceptors"] #MW is the molecular weight, LogP is the solubility, NumHDonors is the hydrogen bond donors, and NumHAcceptors is the hydrogen bond acceptors.
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)
    
    return descriptors

lipinskiDF = lipinski(df.canonical_smiles) #create a new dataframe for molecular descriptors collected by the above.
#print(lipinskiDF) #checkpoint for lipinski function

#combine newDF and lipinskiDF to localize all data
df_combined = pd.concat([newDF, lipinskiDF], axis = 1)
#print(df_combined) #checkpoint for combined dataframes

#convert IC50 to pIC50 --> transforming the IC50 to pIC50 makes the distribution more even
def pIC50(input):
    pIC50 = []
    for i in input['standard_value_norm']:
        molar = i * (10**-9) #converts nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', 1)
    return x

print(df_combined.standard_value.describe()) #checkpoint for pIC50

-np.log10((10**-9) * 100000000)
-np.log10((10**-9) * 10000000000)

def normalValue(input):
    normal = []

    for i in input('standard_value_norm'):
        if i > 100000000:
            i = 100000000
        normal.append(i)

    input['standard_value_norm'] = normal
    x = input.drop('standard_value', 1)
    return x

#def normal = normalValue(df_combined):









#part 1: target protein serach --> collect data set --> bioactivity data --> data processing --> focused data
#part 2: exploratory data analysis
#part 3: Descriptor calculation
#part 4: Modeling 
#part 5: Model Comparison
#part 6: performance comparison 
#part 7: deploy model
#
## Q-star

#focused on coronavirus, specifically MERS-CoV
import pandas as pd
import numpy as np
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski


#Part 1: Target Search --> collecting data sets from ChEMBL database
target = new_client.target
target_query = target.search('coronavirus')
targets = pd.DataFrame.from_dict(target_query)
#print(targets) #checkpoint

selected_target = targets.target_chembl_id[7] #selects for MERS
#print(selected_target) #checkpoint 

#retrieve bioactivity data as the reported IC50 values
bioactivity = new_client.activity
res = bioactivity.filter(target_chembl_id = selected_target).filter(standard_type = "IC50") #double filter function. 1) Searches the database for the selected target then 2) filters for the IC50 values within the database for our selected target.
df = pd.DataFrame.from_dict(res)
df.head(5) #checkpoint --> takes the first 5 / 100 --> add print() to do a full check
df.to_csv('Bioactivity_data.csv', index = False) # save the resulting bioactivity data as a CSV file

df2 = df[df.standard_value.notna()] #removing compounds with missing standard_values
#print(df2) #checkpoint

#pre-processing of bioactivity data
bioactivity_class = []
for i in df2.standard_value:
    if float(i) >=25000:
        bioactivity_class.append("inactive")
    elif float(i) <= 1000:
        bioactivity_class.append("active")
    else:
        bioactivity_class.append('intermediate')

#create a new dataframe that contains unique value for a compound
selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
df3 = df2[selection]
#print(df3) #checkpoint for df3 

#create dataframe to CSV file
df3.to_csv('Bioactivity_preprocessed_data.csv', index = False)

# ---- successfully downloaded chembl database for MERS-CoV ----- #

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









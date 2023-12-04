#part 1: target protein search --> collect data set --> bioactivity data --> data processing --> focused data

#focused on coronavirus, specifically MERS-CoV
import pandas as pd
from chembl_webresource_client.new_client import new_client

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

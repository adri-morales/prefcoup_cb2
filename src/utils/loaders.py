import pandas as pd
import re

def load_signaling_df():

    # Read signling data
    path = 'data/raw/alanin_scanning-cuantitative_values-emax_ec50.csv'
    signaling_df = pd.read_csv(path)

    # Remove mutations that are not alanin scan
    alanin_scan_mask = signaling_df.Mutant.apply(lambda x: bool(re.match('([A-Z])(\d+)$', str(x))))
    signaling_df = signaling_df.loc[alanin_scan_mask]

    # Remove non usefull columns
    cols_of_interst = ['Position','Emax_avg_corr_Gi2_HU210', 'Emax_avg_corr_bArr1_HU210']
    signaling_df = signaling_df[cols_of_interst]

    # Change column names
    signaling_df.columns = ['position', 'gi_emax', 'barr_emax']

    # Set mutant number as index
    signaling_df = signaling_df.set_index('position')

    # Set boolean masks for biased and wt
    biased_mask = (signaling_df.barr_emax == 0) & (signaling_df.gi_emax > 0)
    wt_mask = (signaling_df.barr_emax > 0) & (signaling_df.gi_emax > 0)
        
    # Create a field wih categorical signal information
    signaling_df['profile'] = 2
    signaling_df.loc[wt_mask, 'profile'] = 0
    signaling_df.loc[biased_mask, 'profile'] = 1

    # Keep only the categorical labels
    signaling_df = signaling_df[['profile']]

    # Remove mutants with other signaling profiles
    signaling_df = signaling_df[signaling_df.profile != 2]

    return signaling_df

def load_contacts(path):
    
    # Read tsv file
    contacts = pd.read_csv(path,
                           sep="\t",
                           names=["frame", "type", "atom1",
                                  "atom2", "atom3", "atom4"],
                           skiprows=2,
                           # Ignore atom3 and 4 columns that are only used for water bridges
                           usecols=["frame", "type", "atom1", "atom2"])

    # Remove atom identifier and use only resid and the name of the atom
    for i in [1, 2]:
        # Parse resid and atom name
        atom_data = contacts[f"atom{i}"].str.split(":", expand=True).iloc[:, [
            1, 2, 3]].rename(columns={1: f'resname{i}', 2: f"resid{i}", 3: f"name{i}"})
        # Remove atom column 
        contacts = contacts.drop([f"atom{i}"], axis="columns")
        # Concatenate new data into dataframe
        contacts = pd.concat([contacts, atom_data], axis=1)
        
    return contacts
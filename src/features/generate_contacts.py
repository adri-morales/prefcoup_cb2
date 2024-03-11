""""""

# Standard imports
import argparse
import re
import os
import subprocess
from glob import glob
# Third party imports
import itertools

def get_files(input, traj_format='dcd', top_format='psf'):
    """searches for topology and trajectory files in the input folder and generates a generator 
    of all combinations """
    
    # Find all topology files
    top_files = glob(f"{input}/*.{top_format}")
    # find all trajectory files
    traj_files = glob(f"{input}/*.{traj_format}")
    
    # Generate all combinations of topology and trajectory files
    top_traj_generator = itertools.product(top_files, traj_files)
    
    return top_traj_generator
                  
def compute_contacts(top, traj, sel1='protein', sel2='protein', lib=None, results=None):
    
    """Take pahts names for a topology and trajectory files and generate contact and 
    contact frequency information
    
    This function uses the modules of getcontacts get_dynamic_contacts and get contact_frequencies
    
    parameters:
        top (str): path to a topology file
        traj (str): path to a correponding trajectory file
        lib (str): path to the root of getcontacts. Defaults to src/libraries/getcontacts/
        results (str): path where the results folder will be created. 
    returns:
        None
    """
    # Set default getcontacts path
    if lib == None:
        lib = "src/external/getcontacts/"
    
    # Set the default results path
    if results == None:
        results = "tmp/"
    
    # Extract the name of the system from traj path
    system = traj.split('/')[-2]
    replica = re.match(r'.*(\d)\.', traj.split('/')[-1]).group(1)
    output_folder = f"{results}/{system}"

    
    print(f"Processing system {system}")
    
    # Create the output folder if it does not exist
    if not os.path.exists(output_folder):
        print(f"Creating output folder: {output_folder}")
        os.mkdir(output_folder)
    else:
        print("Output folder already exists")
    
    # Compute contacts
    command_list = [f"{lib}get_dynamic_contacts.py",
                    "--topology",  top,
                    "--trajectory", traj, 
                    "--sele", sel1,
                    "--sele2", sel2,
                    "--output", f"{output_folder}/{replica}contacts.tsv",
                    "--itypes","all",
                    "--cores", "12"]    
    
    print("Computing contacts")    
    subprocess.run(command_list, check=True)
    
    # compute contact frequencies
    print("Computing contact frequencies")
    command = [f"{lib}get_contact_frequencies.py",
                    "--input_files", f"{output_folder}/{replica}contacts.tsv",
                    "--output_file", f"{output_folder}/{replica}frequencies.tsv"]
    subprocess.run(command, check=True)

def main(input, output):
    
    top_traj_generator = get_files(input)
    
    sel1 = 'protein or resname CAN THC'
    sel2 = 'protein or resname CAN THC'
    
    for top, traj, in top_traj_generator:
    
        compute_contacts(top, traj, lib='src/external/getcontacts/', results=output,
                         sel1=sel1, sel2=sel2)

if __name__ == "__main__":
    
    # Parse arguments
    parser = argparse.ArgumentParser()
    
    # add folder argument
    parser.add_argument("--input", type=str, default="data/raw",
                        help="Folder where the data is stored")
    
    # add output argument
    parser.add_argument("--output", type=str, default="data/processed/contacts.csv",
                        help="Path to the output file")
    
    args = dict(parser.parse_args())
    
    main(**args)
    
    # main(input='/gpcr/users/shared/for_adri_veprintsev/canabinol',
    #      output='data/interim/contacts')
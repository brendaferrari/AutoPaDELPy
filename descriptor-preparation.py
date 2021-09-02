import pandas as pd
import sys
import cmd

#Input name and test
if len(sys.argv) <= 1:
    print("One parameter is missing. Please add input file name at following manner: 'python descriptor_preparation.py filename'")
    sys.exit()

filename = sys.argv[1]

#To read molecule dataset (format = smiles, molecule name)
mol_dir = pd.read_csv(filename, names=(['smiles','molecule_id']))

#To save the dataset to use at the descriptors calculation
mol_dir.to_csv('dataset.smi', sep='\t', index=False, header=False)
print(mol_dir)

#To calculate descriptors
## Fingerprint descriptors

from functions.fingerprint_functions import FingerprintFunctions

FingerprintFunctions().cmdloop()

## Bidimensional descriptors

from functions.bidimensional_functions import BidimensionalFunctions

BidimensionalFunctions().cmdloop()

## Tridimensional descriptors

from functions.tridimensional_functions import TridimensionalFunctions

TridimensionalFunctions().cmdloop()

#To move files

from functions.data_processing import moveOutputFile

moveOutputFile()

#To generate descriptors table

from functions.data_processing import saveFile

saveFile()
import cmd
from padelpy import padeldescriptor
import glob

class FingerprintFunctions(cmd.Cmd):
    file = None

    def __init__(self):
        cmd.Cmd.__init__(self)
        self.prompt = '(fingerprint descriptor) '
        self.intro = 'Which fingeprint descriptor do you wish to calculate? Type help or ? to list commands. Write finish if you want to move on to the next descriptor. \n'
        self.completekey='tab'

    def do_AtomPairs2DCount(self, mol_dir):
        """AtomPairs2DCount [mol_dir]
        Calculate the AtomPairs2DFingerprintCount fingerprint"""

        xml_files = glob.glob("functions/descriptors/fingerprint_descriptors/*.xml")
        xml_files.sort()

        FP_list = ['AtomPairs2DCount',
        'AtomPairs2D',
        'EState',
        'CDKextended',
        'CDK',
        'CDKgraphonly',
        'KlekotaRothCount',
        'KlekotaRoth',
        'MACCS',
        'PubChem',
        'SubstructureCount',
        'Substructure']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'AtomPairs2DCount'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        padeldescriptor(mol_dir='dataset.smi',
                        d_file=fingerprint_output_file,
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True)

        return False

    def do_AtomPairs2D(self, mol_dir):
        """AtomPairs2D [mol_dir]
        Calculate the AtomPairs2DFingerprinter fingerprint"""

        xml_files = glob.glob("functions/descriptors/fingerprint_descriptors/*.xml")
        xml_files.sort()

        FP_list = ['AtomPairs2DCount',
        'AtomPairs2D',
        'EState',
        'CDKextended',
        'CDK',
        'CDKgraphonly',
        'KlekotaRothCount',
        'KlekotaRoth',
        'MACCS',
        'PubChem',
        'SubstructureCount',
        'Substructure']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'AtomPairs2D'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        padeldescriptor(mol_dir='dataset.smi',
                        d_file=fingerprint_output_file,
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True)


        return False

    def do_EState(self, mol_dir):
        """EState [mol_dir]
        Calculate the EStateFingerprinter fingerprint"""

        xml_files = glob.glob("functions/descriptors/fingerprint_descriptors/*.xml")
        xml_files.sort()

        FP_list = ['AtomPairs2DCount',
        'AtomPairs2D',
        'EState',
        'CDKextended',
        'CDK',
        'CDKgraphonly',
        'KlekotaRothCount',
        'KlekotaRoth',
        'MACCS',
        'PubChem',
        'SubstructureCount',
        'Substructure']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'EState'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        padeldescriptor(mol_dir='dataset.smi',
                        d_file=fingerprint_output_file,
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True)


        return False

    def do_CDKextended(self, mol_dir):
        """CDKextended[mol_dir]
        Calculate the ExtendedFingerprinter fingerprint"""

        xml_files = glob.glob("functions/descriptors/fingerprint_descriptors/*.xml")
        xml_files.sort()

        FP_list = ['AtomPairs2DCount',
        'AtomPairs2D',
        'EState',
        'CDKextended',
        'CDK',
        'CDKgraphonly',
        'KlekotaRothCount',
        'KlekotaRoth',
        'MACCS',
        'PubChem',
        'SubstructureCount',
        'Substructure']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'CDKextended'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        padeldescriptor(mol_dir='dataset.smi',
                        d_file=fingerprint_output_file,
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True)


        return False

    def do_CDK(self, mol_dir):
        """CDK [mol_dir]
        Calculate the Fingerprinter fingerprint"""

        xml_files = glob.glob("functions/descriptors/fingerprint_descriptors/*.xml")
        xml_files.sort()

        FP_list = ['AtomPairs2DCount',
        'AtomPairs2D',
        'EState',
        'CDKextended',
        'CDK',
        'CDKgraphonly',
        'KlekotaRothCount',
        'KlekotaRoth',
        'MACCS',
        'PubChem',
        'SubstructureCount',
        'Substructure']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'CDK'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        padeldescriptor(mol_dir='dataset.smi',
                        d_file=fingerprint_output_file,
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True)

        return False

    def do_CDKgraphonly(self, mol_dir):
        """CDKgraphonly [mol_dir]
        Calculate the GraphOnlyFingerprinter fingerprint"""

        xml_files = glob.glob("functions/descriptors/fingerprint_descriptors/*.xml")
        xml_files.sort()

        FP_list = ['AtomPairs2DCount',
        'AtomPairs2D',
        'EState',
        'CDKextended',
        'CDK',
        'CDKgraphonly',
        'KlekotaRothCount',
        'KlekotaRoth',
        'MACCS',
        'PubChem',
        'SubstructureCount',
        'Substructure']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'CDKgraphonly'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        padeldescriptor(mol_dir='dataset.smi',
                        d_file=fingerprint_output_file,
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True)


        return False

    def do_KlekotaRothCount(self, mol_dir):
        """KlekotaRothCount [mol_dir]
        Calculate the KlekotaRothFingerprintCount fingerprint"""

        xml_files = glob.glob("functions/descriptors/fingerprint_descriptors/*.xml")
        xml_files.sort()

        FP_list = ['AtomPairs2DCount',
        'AtomPairs2D',
        'EState',
        'CDKextended',
        'CDK',
        'CDKgraphonly',
        'KlekotaRothCount',
        'KlekotaRoth',
        'MACCS',
        'PubChem',
        'SubstructureCount',
        'Substructure']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'KlekotaRothCount'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        padeldescriptor(mol_dir='dataset.smi',
                        d_file=fingerprint_output_file,
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True)


        return False

    def do_KlekotaRoth(self, mol_dir):
        """KlekotaRoth [mol_dir]
        Calculate the KlekotaRothFingerprinter fingerprint"""

        xml_files = glob.glob("functions/descriptors/fingerprint_descriptors/*.xml")
        xml_files.sort()

        FP_list = ['AtomPairs2DCount',
        'AtomPairs2D',
        'EState',
        'CDKextended',
        'CDK',
        'CDKgraphonly',
        'KlekotaRothCount',
        'KlekotaRoth',
        'MACCS',
        'PubChem',
        'SubstructureCount',
        'Substructure']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'KlekotaRoth'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        padeldescriptor(mol_dir='dataset.smi',
                        d_file=fingerprint_output_file,
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True)


        return False

    def do_MACCS(self, mol_dir):
        """MACCS [mol_dir]
        Calculate the MACCSFingerprinter fingerprint"""

        xml_files = glob.glob("functions/descriptors/fingerprint_descriptors/*.xml")
        xml_files.sort()

        FP_list = ['AtomPairs2DCount',
        'AtomPairs2D',
        'EState',
        'CDKextended',
        'CDK',
        'CDKgraphonly',
        'KlekotaRothCount',
        'KlekotaRoth',
        'MACCS',
        'PubChem',
        'SubstructureCount',
        'Substructure']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'MACCS'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        padeldescriptor(mol_dir='dataset.smi',
                        d_file=fingerprint_output_file,
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True)


        return False

    def do_PubChem(self, mol_dir):
        """PubChem [mol_dir]
        Calculate the PubChemFingerprinter fingerprint"""

        xml_files = glob.glob("functions/descriptors/fingerprint_descriptors/*.xml")
        xml_files.sort()

        FP_list = ['AtomPairs2DCount',
        'AtomPairs2D',
        'EState',
        'CDKextended',
        'CDK',
        'CDKgraphonly',
        'KlekotaRothCount',
        'KlekotaRoth',
        'MACCS',
        'PubChem',
        'SubstructureCount',
        'Substructure']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'PubChem'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        padeldescriptor(mol_dir='dataset.smi',
                        d_file=fingerprint_output_file,
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True)


        return False

    def do_SubstructureCount(self, mol_dir):
        """SubstructureCount [mol_dir]
        Calculate the SubstructureFingerprintCount fingerprint"""

        xml_files = glob.glob("functions/descriptors/fingerprint_descriptors/*.xml")
        xml_files.sort()

        FP_list = ['AtomPairs2DCount',
        'AtomPairs2D',
        'EState',
        'CDKextended',
        'CDK',
        'CDKgraphonly',
        'KlekotaRothCount',
        'KlekotaRoth',
        'MACCS',
        'PubChem',
        'SubstructureCount',
        'Substructure']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'SubstructureCount'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        padeldescriptor(mol_dir='dataset.smi',
                        d_file=fingerprint_output_file,
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True)


        return False

    def do_Substructure(self, mol_dir):
        """Substructure [mol_dir]
        Calculate the SubstructureFingerprinter fingerprint"""

        xml_files = glob.glob("functions/descriptors/fingerprint_descriptors/*.xml")
        xml_files.sort()

        FP_list = ['AtomPairs2DCount',
        'AtomPairs2D',
        'EState',
        'CDKextended',
        'CDK',
        'CDKgraphonly',
        'KlekotaRothCount',
        'KlekotaRoth',
        'MACCS',
        'PubChem',
        'SubstructureCount',
        'Substructure']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'Substructure'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        padeldescriptor(mol_dir='dataset.smi',
                        d_file=fingerprint_output_file,
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True)


        return False

    def close(self):
        if self.file:
            self.file.close()
            self.file = None

    def do_finish(self, arg):
        """finish [self, arg]
        Finish running Fingerprint Descriptors and move on to the next descriptor."""

        print('Moving on to the next descriptor...')
        self.close()
        return True

if __name__ == '__main__':
    FingerprintFunctions().cmdloop()
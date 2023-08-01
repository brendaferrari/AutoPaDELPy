import cmd
from padelpy import padeldescriptor
import glob

class BidimensionalFunctions(cmd.Cmd):
    file = None

    def __init__(self):
        cmd.Cmd.__init__(self)
        self.prompt = '(bidimensional descriptor) '
        self.intro = 'Which bidimensional descriptor do you wish to calculate? Type help or ? to list commands. Write finish if you want to move on to the next descriptor.\n'
        self.completekey='tab'
    
    def do_AcidicGroupCount(self, mol_dir):
        """AcidicGroupCount [mol_dir]
        Calculate the AcidicGroupCount fingerprint"""

        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'AcidicGroupCount'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_ALOGP(self, mol_dir):
        """ALOGP [mol_dir]
        Calculate the ALOGP fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'ALOGP'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_AminoAcidCount(self, mol_dir):
        """AminoAcidCount [mol_dir]
        Calculate the AminoAcidCount fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'AminoAcidCount'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_APol(self, mol_dir):
        """APol [mol_dir]
        Calculate the APol fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'APol'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_AromaticAtomsCount(self, mol_dir):
        """AromaticAtomsCount [mol_dir]
        Calculate the AromaticAtomsCount fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'AromaticAtomsCount'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_AromaticBondsCount(self, mol_dir):
        """AromaticBondsCount [mol_dir]
        Calculate the AromaticBondsCount fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'AromaticBondsCount'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_AtomCount(self, mol_dir):
        """AtomCount [mol_dir]
        Calculate the AtomCount fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'AtomCount'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_Autocorrelation(self, mol_dir):
        """Autocorrelation [mol_dir]
        Calculate the Autocorrelation fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'Autocorrelation'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_BaryszMatrix(self, mol_dir):
        """BaryszMatrix [mol_dir]
        Calculate the Autocorrelation fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'BaryszMatrix'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_BasicGroupCount(self, mol_dir):
        """BasicGroupCount [mol_dir]
        Calculate the BasicGroupCount fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'BasicGroupCount'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_BCUT(self, mol_dir):
        """BCUT [mol_dir]
        Calculate the BCUT fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'BCUT'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_BondCount(self, mol_dir):
        """BondCount [mol_dir]
        Calculate the BondCount fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'BondCount'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_BPol(self, mol_dir):
        """BPol [mol_dir]
        Calculate the BPol fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'BPol'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_BurdenModifiedEigenvalues(self, mol_dir):
        """BurdenModifiedEigenvalues [mol_dir]
        Calculate the BurdenModifiedEigenvalues fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'BurdenModifiedEigenvalues'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_CarbonTypes(self, mol_dir):
        """CarbonTypes [mol_dir]
        Calculate the CarbonTypes fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'CarbonTypes'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_ChiChain(self, mol_dir):
        """ChiChain [mol_dir]
        Calculate the ChiChain fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'ChiChain'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_ChiCluster(self, mol_dir):
        """ChiCluster [mol_dir]
        Calculate the ChiCluster fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'ChiCluster'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_ChiPathCluster(self, mol_dir):
        """ChiPathCluster [mol_dir]
        Calculate the ChiPathCluster fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'ChiPathCluster'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_ChiPath(self, mol_dir):
        """ChiPath [mol_dir]
        Calculate the ChiPath fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'ChiPath'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_Constitutional(self, mol_dir):
        """Constitutional [mol_dir]
        Calculate the Constitutional fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'Constitutional'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_Crippen(self, mol_dir):
        """Crippen [mol_dir]
        Calculate the Crippen fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'Crippen'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_DetourMatrix(self, mol_dir):
        """DetourMatrix [mol_dir]
        Calculate the DetourMatrix fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'DetourMatrix'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_EccentricConnectivityIndex(self, mol_dir):
        """EccentricConnectivityIndex [mol_dir]
        Calculate the EccentricConnectivityIndex fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'EccentricConnectivityIndex'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_EStateAtomType(self, mol_dir):
        """EStateAtomType [mol_dir]
        Calculate the EStateAtomType fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'EStateAtomType'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_ExtendedTopochemicalAtom(self, mol_dir):
        """ExtendedTopochemicalAtom [mol_dir]
        Calculate ExtendedTopochemicalAtom fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'ExtendedTopochemicalAtom'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_FMF(self, mol_dir):
        """FMF [mol_dir]
        Calculate FMF fingerprint"""




        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'FMF'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_FragmentComplexity(self, mol_dir):
        """FragmentComplexity [mol_dir]
        Calculate FragmentComplexity fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'FragmentComplexity'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_HBondAcceptorCount(self, mol_dir):
        """HBondAcceptorCount [mol_dir]
        Calculate HBondAcceptorCount fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'HBondAcceptorCount'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_HBondDonorCount(self, mol_dir):
        """HBondDonorCount [mol_dir]
        Calculate HBondDonorCount fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'HBondDonorCount'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_HybridizationRatio(self, mol_dir):
        """HybridizationRatio [mol_dir]
        Calculate HybridizationRatio fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'HybridizationRatio'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_InformationContent(self, mol_dir):
        """InformationContent [mol_dir]
        Calculate InformationContent fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'InformationContent'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_IPMolecularLearning(self, mol_dir):
        """IPMolecularLearning [mol_dir]
        Calculate IPMolecularLearning fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'IPMolecularLearning'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_KappaShapeIndices(self, mol_dir):
        """KappaShapeIndices [mol_dir]
        Calculate KappaShapeIndices fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'KappaShapeIndices'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_KierHallSmarts(self, mol_dir):
        """KierHallSmarts [mol_dir]
        Calculate KierHallSmarts fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'KierHallSmarts'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_LargestChain(self, mol_dir):
        """LargestChain [mol_dir]
        Calculate LargestChain fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'LargestChain'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_LargestPiSystem(self, mol_dir):
        """LargestPiSystem [mol_dir]
        Calculate LargestPiSystem fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'LargestPiSystem'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_LongestAliphaticChain(self, mol_dir):
        """LongestAliphaticChain [mol_dir]
        Calculate LongestAliphaticChain fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'LongestAliphaticChain'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_MannholdLogP(self, mol_dir):
        """MannholdLogP [mol_dir]
        Calculate MannholdLogP fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'MannholdLogP'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_McGowanVolume(self, mol_dir):
        """McGowanVolume [mol_dir]
        Calculate McGowanVolume fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'McGowanVolume'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_MDE(self, mol_dir):
        """MDE [mol_dir]
        Calculate MDE fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'MDE'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_MLFER(self, mol_dir):
        """MLFER [mol_dir]
        Calculate MLFER fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'MLFER'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_PathCount(self, mol_dir):
        """PathCount [mol_dir]
        Calculate PathCount fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'PathCount'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_PetitjeanNumber(self, mol_dir):
        """PetitjeanNumber [mol_dir]
        Calculate PetitjeanNumber fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'PetitjeanNumber'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_RingCount(self, mol_dir):
        """RingCount [mol_dir]
        Calculate RingCount fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'RingCount'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_RotatableBondsCount(self, mol_dir):
        """RotatableBondsCount [mol_dir]
        Calculate RotatableBondsCount fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'RotatableBondsCount'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_RuleOfFive(self, mol_dir):
        """RuleOfFive [mol_dir]
        Calculate RuleOfFive fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'RuleOfFive'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_Topological(self, mol_dir):
        """Topological [mol_dir]
        Calculate Topological fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'Topological'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_TopologicalCharge(self, mol_dir):
        """TopologicalCharge [mol_dir]
        Calculate TopologicalCharge fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'TopologicalCharge'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_TopologicalDistanceMatrix(self, mol_dir):
        """TopologicalDistanceMatrix [mol_dir]
        Calculate TopologicalDistanceMatrix fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'TopologicalDistanceMatrix'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_TPSA(self, mol_dir):
        """TPSA [mol_dir]
        Calculate TPSA fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'TPSA'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_VABC(self, mol_dir):
        """VABC[mol_dir]
        Calculate VABC fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'VABC'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_VAdjMa(self, mol_dir):
        """VAdjMa[mol_dir]
        Calculate VAdjMa fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'VAdjMa'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_WalkCount(self, mol_dir):
        """WalkCount [mol_dir]
        Calculate WalkCount fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'WalkCount'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_Weight(self, mol_dir):
        """Weight [mol_dir]
        Calculate Weight fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'Weight'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_WeightedPath(self, mol_dir):
        """WeightedPath [mol_dir]
        Calculate WeightedPath fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'WeightedPath'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_WienerNumbers(self, mol_dir):
        """WienerNumbers [mol_dir]
        Calculate WienerNumbers fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'WienerNumbers'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_XLogP(self, mol_dir):
        """XLogP [mol_dir]
        Calculate XLogP fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'XLogP'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file, #'Substructure.csv'
                        #descriptortypes='SubstructureFingerprint.xml', 
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_ZagrebIndex(self, mol_dir):
        """ZagrebIndex [mol_dir]
        Calculate ZagrebIndex fingerprint"""



        xml_files = glob.glob("functions/descriptors/bidimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['AcidicGroupCount',
        'ALOGP',
        'AminoAcidCount',
        'APol',
        'AromaticAtomsCount',
        'AromaticBondsCount',
        'AtomCount',
        'Autocorrelation',
        'BaryszMatrix',
        'BasicGroupCount',
        'BCUT',
        'BondCount',
        'BPol',
        'BurdenModifiedEigenvalues',
        'CarbonTypes',
        'ChiChain',
        'ChiCluster',
        'ChiPathCluster',
        'ChiPath',
        'Constitutional',
        'Crippen',
        'DetourMatrix',
        'EccentricConnectivityIndex',
        'EStateAtomType',
        'ExtendedTopochemicalAtom',
        'FMF',
        'FragmentComplexity',
        'HBondAcceptorCount',
        'HBondDonorCount',
        'HybridizationRatio',
        'InformationContent',
        'IPMolecularLearning',
        'KappaShapeIndices',
        'KierHallSmarts',
        'LargestChain',
        'LargestPiSystem',
        'LongestAliphaticChain',
        'MannholdLogP',
        'McGowanVolume',
        'MDE',
        'MLFER',
        'PathCount',
        'PetitjeanNumber',
        'RingCount',
        'RotatableBondsCount',
        'RuleOfFive',
        'Topological',
        'TopologicalCharge',
        'TopologicalDistanceMatrix',
        'TPSA',
        'VABC',
        'VAdjMa',
        'WalkCount',
        'Weight',
        'WeightedPath',
        'WienerNumbers',
        'XLogP',
        'ZagrebIndex']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'ZagrebIndex'

        fingerprint_output_file = ''.join([fingerprint,'.csv'])
        fingerprint_descriptortypes = fp[fingerprint]

        descriptor = padeldescriptor(mol_dir='dataset.smi', 
                        d_file=fingerprint_output_file,
                        descriptortypes= fingerprint_descriptortypes,
                        detectaromaticity=True,
                        standardizenitro=True,
                        standardizetautomers=True,
                        threads=2,
                        removesalt=True,
                        log=True,
                        fingerprints=True,
                        d_2d=True) 

        return False

    def do_finish(self, arg):
        """finish [self, arg]
        Finish running Bidimensional Descriptors and move on to the next descriptor."""

        print('Moving on to the next descriptor...')
        self.close()
        return True

    def close(self):
        if self.file:
            self.file.close()
            self.file = None

if __name__ == '__main__':
    BidimensionalFunctions().cmdloop()
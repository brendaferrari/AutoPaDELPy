import cmd

class TridimensionalFunctions(cmd.Cmd):
    file = None
    
    def __init__(self):
        cmd.Cmd.__init__(self)
        self.prompt = '(tridimensional descriptor) '
        self.intro = 'Which tridimensional descriptor do you wish to calculate? Type help or ? to list commands. Write finish if you want to move on to the preprocessing step.\n'
        self.completekey='tab'

    def do_Autocorrelation3D(self, mol_dir):
        """Autocorrelation3D [mol_dir]
        Calculate the Autocorrelation3D fingerprint"""

        from padelpy import padeldescriptor
        import glob

        xml_files = glob.glob("functions/descriptors/tridimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['Autocorrelation3D',
        'CPSA',
        'GravitationalIndex',
        'LengthOverBreadth',
        'MomentOfInertia',
        'PetitjeanShapeIndex',
        'RDF',
        'WHIM']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'Autocorrelation3D'

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
                        convert3d=True,
                        d_3d=True) 


        return False

    def do_CPSA(self, mol_dir):
        """CPSA [mol_dir]
        Calculate the CPSA fingerprint"""

        from padelpy import padeldescriptor
        import glob

        xml_files = glob.glob("functions/descriptors/tridimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['Autocorrelation3D',
        'CPSA',
        'GravitationalIndex',
        'LengthOverBreadth',
        'MomentOfInertia',
        'PetitjeanShapeIndex',
        'RDF',
        'WHIM']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'CPSA'

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
                        convert3d=True,
                        d_3d=True) 


        return False

    def do_GravitationalIndex(self, mol_dir):
        """GravitationalIndex[mol_dir]
        Calculate the GravitationalIndex fingerprint"""

        from padelpy import padeldescriptor
        import glob

        xml_files = glob.glob("functions/descriptors/tridimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['Autocorrelation3D',
        'CPSA',
        'GravitationalIndex',
        'LengthOverBreadth',
        'MomentOfInertia',
        'PetitjeanShapeIndex',
        'RDF',
        'WHIM']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'GravitationalIndex'

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
                        convert3d=True,
                        d_3d=True) 


        return False

    def do_LengthOverBreadth(self, mol_dir):
        """LengthOverBreadth [mol_dir]
        Calculate the LengthOverBreadth fingerprint"""

        from padelpy import padeldescriptor
        import glob

        xml_files = glob.glob("functions/descriptors/tridimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['Autocorrelation3D',
        'CPSA',
        'GravitationalIndex',
        'LengthOverBreadth',
        'MomentOfInertia',
        'PetitjeanShapeIndex',
        'RDF',
        'WHIM']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'LengthOverBreadth'

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
                        convert3d=True,
                        d_3d=True) 


        return False

    def do_MomentOfInertia(self, mol_dir):
        """MomentOfInertia [mol_dir]
        Calculate the MomentOfInertia fingerprint"""

        from padelpy import padeldescriptor
        import glob

        xml_files = glob.glob("functions/descriptors/tridimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['Autocorrelation3D',
        'CPSA',
        'GravitationalIndex',
        'LengthOverBreadth',
        'MomentOfInertia',
        'PetitjeanShapeIndex',
        'RDF',
        'WHIM']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'MomentOfInertia'

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
                        convert3d=True,
                        d_3d=True) 


        return False

    def do_PetitjeanShapeIndex(self, mol_dir):
        """PetitjeanShapeIndex [mol_dir]
        Calculate the PetitjeanShapeIndex fingerprint"""

        from padelpy import padeldescriptor
        import glob

        xml_files = glob.glob("functions/descriptors/tridimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['Autocorrelation3D',
        'CPSA',
        'GravitationalIndex',
        'LengthOverBreadth',
        'MomentOfInertia',
        'PetitjeanShapeIndex',
        'RDF',
        'WHIM']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'PetitjeanShapeIndex'

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
                        convert3d=True,
                        d_3d=True) 


        return False

    def do_RDF(self, mol_dir):
        """RDF [mol_dir]
        Calculate the RDF fingerprint"""

        from padelpy import padeldescriptor
        import glob

        xml_files = glob.glob("functions/descriptors/tridimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['Autocorrelation3D',
        'CPSA',
        'GravitationalIndex',
        'LengthOverBreadth',
        'MomentOfInertia',
        'PetitjeanShapeIndex',
        'RDF',
        'WHIM']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'RDF'

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
                        convert3d=True,
                        d_3d=True) 


        return False

    def do_WHIM(self, mol_dir):
        """WHIM [mol_dir]
        Calculate the WHIM fingerprint"""

        from padelpy import padeldescriptor
        import glob

        xml_files = glob.glob("functions/descriptors/tridimensional_descriptors/*.xml")
        xml_files.sort()    

        FP_list = ['Autocorrelation3D',
        'CPSA',
        'GravitationalIndex',
        'LengthOverBreadth',
        'MomentOfInertia',
        'PetitjeanShapeIndex',
        'RDF',
        'WHIM']

        fp = dict(zip(FP_list, xml_files))

        fingerprint = 'WHIM'

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
                        convert3d=True,
                        d_3d=True) 


        return False

    def do_finish(self, arg):
        """finish [self, arg]
        Finish running Tridimensional Descriptors and move on to the next descriptor."""

        print('Moving on to data preprocessing...')
        self.close()
        return True

    def close(self):
        if self.file:
            self.file.close()
            self.file = None

if __name__ == '__main__':
    TridimensionalFunctions().cmdloop()
Amino Acid Fragmenter 

This script aims to fragment protein systems for the purposes of F-ISAPT Analysis. 

input: well-formed pdb

output: Psi4 input file, fA.dat and fB.dat fragmentation files

Instructions:

    To use (in PyMOL):
    1)  Open pdb in PyMOL
    2)  Run aa_chopper.py in PyMOL (Note: this step can also be achieved by opening the script and .pdb 
        simultaneously, perhaps through aliasing)
    3)  Type 'chop'
    
    This will produce:
    1)  A colored fragmentation scheme along with fragment names
    2)  A working input file with some sensible defaults
    3)  fsapt/fA.dat and fsapt/fB.dat files for post-analysis
            -Warning: it will overwrite any fA.dat or fB.dat files living in fsapt/ in the CWD
            -PDB file must have '.pdb' extension
            -Input file produced will have the same name as pdb, with extension '.in'
            -Assumes there is only one ligand, and it is its own fragment
            -Assumes solvent does not carry charge and belongs to monomer C
            -Assumes well-formed pdb supplied with correct atom types for atoms participating in peptide bond
            -Fragment names will be as follows:
                N-terminus caps:    Chain_AminoAcidAbbreviation+ResidueNumber_NTC
                C-terminus caps:    Chain_AminoAcidAbbreviation+ResidueNumber_CTC
                Sidechain caps:     Chain_AminoAcidAbbreviation+ResidueNumber_CS
                Sidechains:         Chain_AminoAcidAbbreviation+ResidueNumber_SC
                Peptide bonds:      Chain_ResidueNumberOfAA1_ResidueNumberOfAA2_PEPT
                Free Amino Acids:   Chain_AminoAcidAbbreviation+ResidueNumber_FREE
                Ligand:             LIG_ResidueName_ResidueNumber
                Solvent:            SOLV_ResidueName_ResidueNumber
    4)  A PyMOL file named the same as the given PDB file, but with .pymol extension
            -This will produce the same visualization as described above when opened with PyMOL
            -This can be done by typing 'pymol -l NAME_OF_FILE.pymol'
    
    To use (in terminal):
    1) Type 'python aa_chopper.py FILENAME.pdb'
    
    This will produce: 
    1)  A working input file with some sensible defaults
    2)  fsapt/fA.dat and fsapt/fB.dat files for post-analysis
            -Named as described above
    3)  A PyMOL file named the same as the given PDB file, but with .pymol extension
            -This will produce the same visualization as described above when opened with PyMOL
            -This can be done by typing 'pymol -l NAME_OF_FILE.pymol'

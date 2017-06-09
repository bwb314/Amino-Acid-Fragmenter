import math, os, sys
from pymol import cmd

# Amino Acid dictionary
aa_dict = {'ALA': 'A' , 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN':
'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

#this function only works if proposed filename has one period and one underscore 
def pick_file_name(proposed_filename, filnames):
    psplit = proposed_filename.split('.')
    base = '_'.join(psplit[0].split('_')[:-1])
    ext = psplit[1]
    count = 0
    for fil in filnames:
        if fil == proposed_filename: 
            count += 1
    filename = base+'_'+str(count)+'.'+ext
    return filename

#color fragments based on fragmentation from bfs
def color_frags(fA, single_color = ''):
    warm = []
    x = math.ceil((len(fA)/8.)**(1./3.))
    for li in range(int(2*x)):
        i = 1.0 - li*(1.0-0.4)/(2*x-1)
        for lj in range(int(4*x)):
            j = lj*(1.0-0.0)/(4*x-1)
            for lk in range(int(x)):
                safe = x
                if x != 1: safe = x - 1
                k = lk*(0.25-0.0)/(safe)
                warm.append([i,j,k])
    
    for frag in range(len(fA)):
        selection = ''
        for atom in range(len(fA[frag][1:])):
            key = str(fA[frag][1+atom])
            selection += 'rank '+key+' '                
        cmd.select("("+fA[frag][0]+")", selection)
        #linspace through color vector
        safefA = len(fA)
        if len(fA) != 1: safefA = len(fA) - 1
        fragind = int(math.floor(frag*((len(warm))-1)/safefA))
        if single_color == '': color = warm[fragind]
        else: color = single_color
        cmd.set_color("A_"+fA[frag][0], color)
        cmd.color("A_"+fA[frag][0],"("+fA[frag][0]+")")



def frag_type(residue):
    is_n = False
    is_c = False
    for atom in residue:
        atom_type = atom[1]
        if atom_type == 'N': is_n = True
        if atom_type == 'C': is_c = True
    if is_n and is_c: return 'SIDECHAIN'
    elif is_n: return 'C_CAP'
    elif is_c: return 'N_CAP'
    else: return 'ISOLATED'

def fA_stuffer(residue,f_type,fA):
    fA_entry = []
    if f_type == 'N_CAP':
        name = residue[0][3]+'_'+residue[0][4]+'_NTC'
    elif f_type == 'C_CAP':
        name = residue[0][3]+'_'+residue[0][4]+'_CTC'
    elif f_type == 'SIDECHAIN':
        name = residue[0][3]+'_'+aa_dict[residue[0][2]]+residue[0][4]+'_SC'
    else:
        name = residue[0][3]+'_'+aa_dict[residue[0][2]]+residue[0][4]+'_FREE'
    fA_entry.append(name)
    for atom in residue: 
        if atom[1] not in ['C','O','N','H']: fA_entry.append(atom[0])
    if len(fA_entry) == 1:
        for atom in residue: fA_entry.append(atom[0])
    fA.append(fA_entry)

def gather_peptides(residue,peptides):
    a_types = ['H','N','C','O']
    pep = []
    for atom in residue: 
        if atom[1] in a_types: pep.append(atom)
    #sort in for stitching
    for atom in a_types:
        for atom_s in pep:
            if atom_s[1] == atom: peptides.append(atom_s)
    
def peptide_stuffer(peptides,fA):
    #num_of_peps = len(peptides)/4
    #for pep_ind in range(num_of_peps):
    #    current_ind = pep_ind * 4
    #    name = peptides[current_ind][4]+'_'+peptides[current_ind+2][4]+'_pept'
    #    fA.append([name,peptides[current_ind][0],peptides[current_ind+1][0],peptides[current_ind+2][0],peptides[current_ind+3][0]])
    pep = []
    # take advantage of sortedness
    for pep_atom in peptides:
        pep.append(pep_atom[0]) 
        if pep_atom[1] == 'C': num1 = pep_atom[4]
        if pep_atom[1] == 'N': 
            num2 = pep_atom[4]
            if num1 == '': name = pep_atom[3]+'_'+num2+'_NTC'
            else: name = pep_atom[3]+'_'+num1+'_'+num2+'_pept'
            pep.insert(0,name)
            fA.append(pep)
            pep = []
            num1 = ''
            num2 = ''
    
def write_frag_file(fA,name):
    fil = open('fsapt/'+name+'.dat','w')
    for entry in fA:
        fil.write(entry[0]+' ')
        for item in entry[1:]: fil.write(str(item+1)+' ')
        fil.write('\n')
    fil.close()

def write_input(protein,ligand, protein_charge, lig_charge, solvent, fil_name):
    inp_fil = open(fil_name.split('.')[0]+'.in','w')
    # psi4 input file sections
    inp_fil.write('molecule {\n'+str(protein_charge)+' 1\n')
    for atom in protein: inp_fil.write(atom[0]+' '+' '.join(atom[1])+'\n')
    inp_fil.write('--\n'+str(lig_charge)+' 1\n')
    for atom in ligand: inp_fil.write(atom[0]+' '+' '.join(atom[1])+'\n')
    if len(solvent) != 0:
        inp_fil.write('--\n0 1\n')
        for atom in solvent: inp_fil.write(atom[0]+' '+' '.join(atom[1])+'\n')
    inp_fil.write('''units angstrom\nno_com\nno_reorient\nsymmetry c1
    }\n\nmemory 30000 MB\n
set {
  basis         jun-cc-pvdz
  df_basis_scf  jun-cc-pvdz-jkfit
  df_basis_sapt jun-cc-pvdz-ri
  scf_type df
  guess sad
  minao_basis cc-pvtz-minao
  maxiter 100
  ints_tolerance 0.0
}\n\nenergy('fisapt0')\n\n''')
    inp_fil.close()

def chop():
    #READ PDB and LIGAND
    fA = []
    fil_name = cmd.get_names("all")[0]+'.pdb'
    pdb_fil = []
    residue = []
    prevres = ''
    ind = 0
    protein_inp = []
    protein_charge = 0
    for lin in open(fil_name,'r'):
        if 'ATOM' not in lin and prevres == '': continue
        try:
            a_type = lin.split()[2]
            res_type = lin.split()[3]
            chain_name = lin.split()[4]
            res_num = lin.split()[5]
            coords = lin.split()[6:9]
            atom = lin.split()[-1]
            if '-' in atom: 
                charge = atom [-2]
                atom = atom[:-2]
                protein_charge -= int(charge)
            if '+' in atom: 
                charge = atom [-2]
                atom = atom[:-2]
                protein_charge += int(charge)
            protein_inp.append([atom,coords])
            if prevres == '': prevres = res_num
        except: 
            pdb_fil.append(residue)
            prevres = ''
            continue
        atom = [ind,a_type,res_type,chain_name,res_num]
        ind += 1
        if res_num != prevres: 
            pdb_fil.append(residue)
            residue = []
            prevres = res_num
        residue.append(atom)
    
    ligands = {}
    ligand_inp = []
    lig_ind = 0
    lig_charge = 0
    solv_inp = []
    #assumes one ligand at the moment
    for lin in open(fil_name,'r'):
        if 'ATOM' in lin: lig_ind += 1 
        if 'HETATM' not in lin: continue
        lignum = lin.split()[5]
        ligchain = lin.split()[4]
        if len(ligchain) > 1: lignum = ligchain[1:]
        ligtype = lin.split()[3]+ '_' + lignum
        coords = lin.split()[6:9]
        atom = lin.split()[-1]
        print ligtype
        if '-' in atom: 
            charge = atom [-2]
            atom = atom[:-2]
            lig_charge -= int(charge)
        if '+' in atom: 
            charge = atom [-2]
            atom = atom[:-2]
            lig_charge += int(charge)
        if ligands == {}:
            ligands[ligtype] = ['LIG_'+ligtype, lig_ind]
            lig_ind += 1
            ligand_inp.append([atom,coords])
            continue
        elif ligtype not in ligands:
            ligands[ligtype] = ['SOLV_'+ligtype, lig_ind]
            lig_ind += 1
            solv_inp.append([atom,coords])
            continue
        if ligtype in ligands:
            ligands[ligtype].append(lig_ind)
            lig_ind += 1
        if 'LIG' in ligands[ligtype][0]: 
            ligand_inp.append([atom,coords])
        else: 
            solv_inp.append([atom,coords])
    
    fB = [[]]
    for i in ligands.keys(): 
        if not ligands[i][0].startswith('LIG'): continue
        fB[0].append(ligands[i][0])
        for j in ligands[i][1:]: fB[0].append(j)
    fC = []
    if len(ligands.keys()) > 1:
        for i in ligands.keys():
            if ligands[i][0].startswith('LIG'): continue
            entry = [ligands[i][0]]
            for j in ligands[i][1:]: entry.append(j)
            fC.append(entry)
    
    peptides = []
    for residue in pdb_fil:
        f_type = frag_type(residue) 
        fA_stuffer(residue,f_type,fA) 
        gather_peptides(residue,peptides)
    
    #print peptides
    peptide_stuffer(peptides,fA)
    
    color_frags(fA)
    color_frags(fB, '[0.5,0.5,0.5]')
    color_frags(fC, '[0.0,0.0,1.0]')
    cmd.show('sticks')
    cmd.label('all','elem')
    try: os.mkdir('fsapt')
    except: print '\n\nfsapt/ already exists. This will likely overwrite old fA.dat and fB.dat.\n\n'
    write_frag_file(fA,'fA')
    write_frag_file(fB,'fB')
    write_input(protein_inp,ligand_inp, protein_charge, lig_charge, solv_inp, fil_name)

cmd.extend("chop",chop)

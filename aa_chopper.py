import math, os, sys
pymol = True
try: from pymol import cmd
except: pymol = False

# Amino Acid dictionary
aa_dict = {'ALA': 'A' , 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN':
'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

#color fragments based on fragmentation from bfs
def color_frags(fA, pymol_file, single_color = ''):
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
        #linspace through color vector
        safefA = len(fA)
        if len(fA) != 1: safefA = len(fA) - 1
        fragind = int(math.floor(frag*((len(warm))-1)/safefA))
        if single_color == '': color = warm[fragind]
        else: color = single_color
        pymol_file.write('cmd.select("('+fA[frag][0]+')", "'+selection+'")\n')
        pymol_file.write('cmd.set_color("A_'+fA[frag][0]+'", '+str(color)+')\n')
        pymol_file.write('cmd.color("A_'+fA[frag][0]+'","('+fA[frag][0]+')")\n')
        if pymol: 
            cmd.select("("+fA[frag][0]+")", selection)
            cmd.set_color("A_"+fA[frag][0], color)
            cmd.color("A_"+fA[frag][0],"("+fA[frag][0]+")")

def write_frag_file(fA,name):
    fil = open('fsapt/'+name+'.dat','w')
    for entry in fA:
        fil.write(entry[0]+' ')
        for item in entry[1:]: fil.write(str(int(item)+1)+' ')
        fil.write('\n')
    fil.close()

def write_input(protein,ligand, protein_charge, lig_charge, solvent, file_name, solv_mon):
    inp_fil = open(file_name.split('.')[0]+'.in','w')
    # psi4 input file sections
    inp_fil.write('molecule {\n'+str(protein_charge)+' 1\n')
    for atom in protein: inp_fil.write(atom[0]+' '+' '.join(atom[1])+'\n')
    if solv_mon.lower() == "a" and len(solvent) != 0:
        for atom in solvent: inp_fil.write(atom[0]+' '+' '.join(atom[1])+'\n')
    inp_fil.write('--\n'+str(lig_charge)+' 1\n')
    for atom in ligand: inp_fil.write(atom[0]+' '+' '.join(atom[1])+'\n')
    if solv_mon.lower() == "b" and len(solvent) != 0:
        for atom in solvent: inp_fil.write(atom[0]+' '+' '.join(atom[1])+'\n')
    if solv_mon.lower() == "c" and len(solvent) != 0:
        inp_fil.write('--\n0 1\n')
        for atom in solvent: inp_fil.write(atom[0]+' '+' '.join(atom[1])+'\n')
    inp_fil.write('''units angstrom\nno_com\nno_reorient\nsymmetry c1
    }\n\nmemory 30000 MB\n
set {
  basis         jun-cc-pvdz
  df_basis_scf  jun-cc-pvdz-jkfit
  df_basis_sapt jun-cc-pvdz-ri
  freeze_core   true
  scf_type df
  guess sad
  minao_basis cc-pvtz-minao
  maxiter 100
  ints_tolerance 0.0
}\n\nenergy('fisapt0')\n\n''')
    inp_fil.close()

def dist(a,b):
    a_coords = [float(x) for x in a['coords'].split()]
    b_coords = [float(x) for x in b['coords'].split()]
    xdist = (a_coords[0]-b_coords[0])**2
    ydist = (a_coords[1]-b_coords[1])**2
    zdist = (a_coords[2]-b_coords[2])**2
    dist = (xdist+ydist+zdist)**0.5
    return dist

cov_rad = {   'H' : 0.37, 'C' : 0.77, 'O' : 0.73, 'N' : 0.75, 'F' : 0.71,
  'P' : 1.10, 'S' : 1.03, 'Cl': 0.99, 'Br': 1.14, 'I' : 1.33, 'He': 0.30,
  'Ne': 0.84, 'Ar': 1.00, 'Li': 1.02, 'Be': 0.27, 'B' : 0.88, 'Na': 1.02,
  'Mg': 0.72, 'Al': 1.30, 'Si': 1.18, 'K' : 1.38, 'Ca': 1.00, 'Sc': 0.75,
  'Ti': 0.86, 'V' : 0.79, 'Cr': 0.73, 'Mn': 0.67, 'Fe': 0.61, 'Co': 0.64,
  'Ni': 0.55, 'Cu': 0.46, 'Zn': 0.60, 'Ga': 1.22, 'Ge': 1.22, 'As': 1.22,
  'Se': 1.17, 'Br': 1.14, 'Kr': 1.03}

def bound(a,b):
    if cov_rad[a['element']]+cov_rad[b['element']] >=  dist(a,b): return True
    else: return False

def intersect(frag,frags):
    intersects = set() 
    for f in frags:
        if frag & f == set([]): continue
        intersects = intersects | frag | f
        frags.remove(f) 
    return intersects

def read_pdb(file_name):
    pdb_file = []
    for line in open(file_name,'r'):
        if 'ATOM' != line[:4] and 'HETATM' != line[:6]: continue
        atom_class = line[:6].split()[0]
        #starting from 0
        atom_index = str(len(pdb_file))
        atom_type = line[12:16].split()[0]
        residue_name = line[17:20].split()[0]
        chain = line[21].split()[0]
        residue_number = line[22:26].split()[0]
        coords = ' '.join(line[30:54].split())
        element = line[76:78].split()[0]
        try: charge = line[78:80].split()[0]
        except: charge = '0'
        atom_info = {'atom_class':atom_class,'atom_index':atom_index,'atom_type':atom_type,'residue_name':residue_name,'chain':chain,'residue_number':residue_number,'coords':coords,'element':element,'charge':charge}
        pdb_file.append(atom_info)
  
    return pdb_file
    
def bfs(peptides, residues):
    translator = {}
    for pep_dict in peptides: translator[str(pep_dict)] = pep_dict

    peptide_frags = {}
    minifrags = []
    #make sets of bonded atoms
    for i in range(len(peptides)):
        for j in range(len(peptides)):
            if j>=i: continue
            if bound(peptides[i],peptides[j]): minifrags.append(set([str(peptides[i]),str(peptides[j])]))
    
    while minifrags != []:
        one = intersect(minifrags[0],minifrags)
        last = []
        while one != set([]):
            one = intersect(one,minifrags)
            if one != set([]): last = one
        if last != []:
            dict_last = [translator[x] for x in last]
            nums = sorted(set([x['residue_number'] for x in dict_last]))
            if len(nums) == 1:
                print 'Detecting that terminating peptide is present.'
                resname = dict_last[0]['chain'] + dict_last[0]['residue_number']
                for a in dict_last: residues[resname].append(a)
                continue
            peptide_name = dict_last[0]['chain']+'_'+nums[0]+'_'+nums[1]+'_PEPT'
            peptide_frags[peptide_name] = dict_last
    
    return peptide_frags

def interpret(pdb_file):
    peptides = []
    residues = {}
    ligand = []
    solvent = {}
    protein_charge = 0
    ligand_charge = 0
    protein_for_inp = []
    ligand_for_inp = []
    solvent_for_inp = []
    for atom in pdb_file:
        if atom['atom_class'] == 'ATOM':
            if atom['atom_type'] in ['C','N','H','O']:
                peptides.append(atom)
            else:
                res_id = atom['chain']+atom['residue_number']
                if res_id in residues:
                    residues[res_id].append(atom)
                else: residues[res_id] = [atom]
            if '-' in atom['charge']:
                protein_charge -= int(atom['charge'].split('-')[0])
            elif '+' in atom['charge']:
                protein_charge += int(atom['charge'].split('+')[0])
            protein_for_inp.append([atom['element'],atom['coords'].split()])
        else:
            if ligand == []:
                ligand_id = atom['residue_number']
                ligand.append(atom)
                ligand_for_inp.append([atom['element'],atom['coords'].split()])
            elif atom['residue_number'] == ligand_id:
                ligand.append(atom)
                ligand_for_inp.append([atom['element'],atom['coords'].split()])
            else:
                if atom['residue_number'] in solvent:
                    solvent[atom['residue_number']].append(atom)
                    solvent_for_inp.append([atom['element'],atom['coords'].split()])
                else: 
                    solvent[atom['residue_number']] = [atom]
                    solvent_for_inp.append([atom['element'],atom['coords'].split()])
            if '-' in atom['charge']:
                ligand_charge -= int(atom['charge'].split('-')[0])
            elif '+' in atom['charge']:
                ligand_charge += int(atom['charge'].split('+')[0])
    
    return residues, peptides, ligand, solvent, protein_charge, ligand_charge, protein_for_inp, ligand_for_inp, solvent_for_inp
    
def fragment(peptide_frags, residues, ligand, solvent, fA, fB, fC, solv_mon):
    segs = []
    seg = {}
    #assumes that peptide bond type atoms form peptide bonds
    for pepind in range(len(sorted(peptide_frags.keys()))):
        pepkey = sorted(peptide_frags.keys())[pepind]
        pepchain = pepkey[0]
        pepnum1 = pepkey.split('_')[1] 
        pepnum2 = pepkey.split('_')[2]
        res1 = pepchain+pepnum1 
        res2 = pepchain+pepnum2 
        if seg == {}:
            try: seg[res1+'_NTC'] = residues[res1]
            except: continue 
            entry = [x['atom_index'] for x in residues[res1]]
            name = res1[0]+'_'+aa_dict[residues[res1][0]['residue_name']]+res1[1:]
            capname = name+'_NTC'
            atoms = [x['element'] for x in residues[res1]]
            #is cap a methyl group? If not capping sidechain
            if atoms.count('C') != 1: capname = name+'_CS'
            if atoms.count('H') != 3: capname = name+'_CS'
            if len(atoms) != 4: capname = name+'_CS'
            entry.insert(0,capname)
            fA.append(entry) 
            del residues[res1]
        else: 
            seg[res1+'_SC'] = residues[res1]
            entry = [x['atom_index'] for x in residues[res1]]
            name = res1[0]+'_'+aa_dict[residues[res1][0]['residue_name']]+res1[1:]
            entry.insert(0,name+'_SC')
            fA.append(entry) 
            del residues[res1] 
        seg[pepkey] = peptide_frags[pepkey]
        entry = [x['atom_index'] for x in list(peptide_frags[pepkey])]
        entry.insert(0,pepkey)
        fA.append(entry) 
        try: pepnum3 = sorted(peptide_frags.keys())[pepind+1].split('_')[2]
        except: pepnum3 = ''
        nextpept = pepchain+'_'+pepnum2+'_'+pepnum3+'_PEPT'
        if nextpept not in peptide_frags.keys():
            try: 
                seg[res2+'_CTC'] = residues[res2]
                entry = [x['atom_index'] for x in residues[res2]]
                name = res2[0]+'_'+aa_dict[residues[res2][0]['residue_name']]+res2[1:]
                capname = name+'_NTC'
                atoms = [x['element'] for x in residues[res2]]
                #is cap a methyl group? If not capping sidechain
                if atoms.count('C') != 1: capname = name+'_CS'
                if atoms.count('H') != 3: capname = name+'_CS'
                if len(atoms) != 4: capname = name+'_CS'
                entry.insert(0,capname)
                fA.append(entry) 
                del residues[res2]
                segs.append(seg)
                seg = {}
            # assume terminal peptide
            except: continue 
   
    for free in sorted(residues.keys()):
        entry = [x['atom_index'] for x in residues[free]]
        name = free[0]+'_'+aa_dict[residues[free][0]['residue_name']]+free[1:]
        entry.insert(0,name+'_FREE')
        fA.append(entry)
        del residues[free]

    #assumes one ligand residue type
    entry = [x['atom_index'] for x in ligand]
    name = 'LIG_'+ligand[0]['residue_name']+'_'+ligand[0]['residue_number']
    entry.insert(0,name)
    fB.append(entry)

    solv_keys = {'a':fA,'b':fB,'c':fC}
    #if different residue type than ligand is solvent
    #assumes solvent goes in c
    for solv_mol in sorted(solvent.keys()):
        entry = [x['atom_index'] for x in solvent[solv_mol]]
        name = 'SOLV_'+solvent[solv_mol][0]['residue_name']+'_'+solvent[solv_mol][0]['residue_number']
        entry.insert(0,name)
        solv_keys[solv_mon.lower()].append(entry) 

def chop(solv_mon = 'a'):
   
    if pymol: 
        file_name = cmd.get_names("all")[0]+'.pdb'
        cmd.show('sticks') 
    else: 
        file_name = ""
        for arg in sys.argv: 
            if arg.endswith('.pdb'): file_name = arg
        if file_name == "": 
            print 'Must provide filename if running from terminal!'
            sys.exit()

    pdb_file = read_pdb(file_name) 

    #assumes solvent does not carry charge
    residues, peptides, ligand, solvent, protein_charge, ligand_charge, protein_for_inp, ligand_for_inp, solvent_for_inp = interpret(pdb_file)   
 
    peptide_frags = bfs(peptides, residues)

    fA, fB, fC = [], [], []
    fragment(peptide_frags, residues, ligand, solvent, fA, fB, fC, solv_mon)

    pymol_file = open(file_name.split('.pdb')[0]+'.pymol','w')
    pymol_file.write("from pymol import cmd\n")
    pymol_file.write("cmd.load('"+file_name+"')\n")
    pymol_file.write("cmd.show('sticks')\n")
    color_frags(fA, pymol_file) 
    color_frags(fB, pymol_file, [0.5,0.5,0.5]) 
    color_frags(fC, pymol_file, [0.0,0.0,1.0]) 
    pymol_file.close()
    try: os.mkdir('fsapt/')
    except: print 'Overwriting old fA.dat and fB.dat'
    write_frag_file(fA,'fA')    
    write_frag_file(fB,'fB')
    write_input(protein_for_inp, ligand_for_inp, protein_charge, ligand_charge, solvent_for_inp, file_name, solv_mon)

if pymol: cmd.extend("chop",chop)
else: 
    args = ' '.join(sys.argv)
    if 'solv_mon' in args and '=' in args:
        solv_choice = args.split('=')[1].split()[0].lower()
    else: solv_choice = 'a'
    chop(solv_mon=solv_choice)

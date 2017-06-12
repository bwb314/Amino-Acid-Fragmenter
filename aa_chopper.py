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



def write_frag_file(fA,name):
    fil = open('fsapt/'+name+'.dat','w')
    for entry in fA:
        fil.write(entry[0]+' ')
        for item in entry[1:]: fil.write(str(int(item)+1)+' ')
        fil.write('\n')
    fil.close()

def write_input(protein,ligand, protein_charge, lig_charge, solvent, file_name):
    inp_fil = open(file_name.split('.')[0]+'.in','w')
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

def dist(a,b):
    a_coords = [float(x) for x in a[6].split()]
    b_coords = [float(x) for x in b[6].split()]
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
    if cov_rad[a[7]]+cov_rad[b[7]] >=  dist(a,b): return True
    else: return False

def intersect(frag,frags):
    intersects = set() 
    for f in frags:
        if frag & f == set([]): continue
        intersects = intersects | frag | f
        frags.remove(f) 
    return intersects

def str2list(string):
    string = string.replace('[','')
    string = string.replace(']','')
    string = string.replace("'",'')
    arr = string.split(',')
    return arr

def chop():
    pdb_file = []
    file_name = cmd.get_names("all")[0]+'.pdb'
    for line in open(file_name,'r'):
        if 'ATOM' not in line and 'HETATM' not in line: continue
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
        atom_info = [atom_class,atom_index,atom_type,residue_name,chain,residue_number,coords,element,charge]
        #atom_info = [atom_class,atom_index,atom_type,residue_name,chain,residue_number,coords,element,charge]
        #                   0         1          2       3            4     5               6    7      8
        pdb_file.append(atom_info)

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
        atom_class = atom[0]
        atom_type = atom[2]
        residue_number = atom[5]
        chain = atom[4]
        charge = atom[8]
        element = atom[7]
        coords = atom[6].split()
        if atom_class == 'ATOM':
            if atom_type in ['C','N','H','O']:
                peptides.append(atom)
            else:
                res_id = chain+residue_number
                if res_id in residues:
                    residues[res_id].append(atom)
                else: residues[res_id] = [atom]
            if '-' in charge:
                protein_charge -= int(charge.split('-')[0])
            elif '+' in charge:
                protein_charge += int(charge.split('+')[0])
            protein_for_inp.append([element,coords])
        else:
            if ligand == []:
                ligand_id = residue_number
                ligand.append(atom)
                ligand_for_inp.append([element,coords])
            elif residue_number == ligand_id:
                ligand.append(atom)
                ligand_for_inp.append([element,coords])
            else:
                if residue_number in solvent:
                    solvent[residue_number].append(atom)
                    solvent_for_inp.append([element,coords])
                else: 
                    solvent[residue_number] = [atom]
                    solvent_for_inp.append([element,coords])
            #assumes only ligand can be charged
            if '-' in charge:
                ligand_charge -= int(charge.split('-')[0])
            elif '+' in charge:
                ligand_charge += int(charge.split('+')[0])

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
        #atom_info = [atom_class,atom_index,atom_type,residue_name,chain,residue_number,coords,element]
        #                   0         1          2       3            4     5               6    7
            list_last = [str2list(x) for x in list(last)]
            nums = sorted(set([x[5] for x in list_last]))
            if len(nums) == 1:
                print 'Terminating peptide'
                resname = list_last[0][4].split()[0] + list_last[0][5].split()[0]
                for a in list_last:
                    new_arr = []
                    for entry in a:
                        if entry == 'ATOM': 
                            new_arr.append(entry)
                        else:
                            new_arr.append(entry[1:])
                    residues[resname].append(new_arr)
                continue
            peptide_name = list_last[0][4].split()[0]+'_'+nums[0].split()[0]+'_'+nums[1].split()[0]+'_PEPT'
            peptide_frags[peptide_name] = last

    #residues = {}
    #ligand = []
    #solvent = {}
    #peptide_frags = {}

    segs = []
    seg = {}
    fA = []
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
            except:
                continue 
                #print residues.keys()
                #print res1,'not present'
                #sys.exit()
            entry = [x[1] for x in residues[res1]]
            name = res1[0]+'_'+aa_dict[residues[res1][0][3]]+res1[1:]
            entry.insert(0,name+'_NTC')
            fA.append(entry) 
            del residues[res1]
        else: 
            seg[res1+'_SC'] = residues[res1]
            entry = [x[1] for x in residues[res1]]
            name = res1[0]+'_'+aa_dict[residues[res1][0][3]]+res1[1:]
            entry.insert(0,name+'_SC')
            fA.append(entry) 
            del residues[res1] 
        seg[pepkey] = peptide_frags[pepkey]
        entry = [str2list(x)[1].split()[0] for x in list(peptide_frags[pepkey])]
        entry.insert(0,pepkey)
        fA.append(entry) 
        #del peptide_frags[pepkey]
        try: pepnum3 = sorted(peptide_frags.keys())[pepind+1].split('_')[2]
        except: pepnum3 = ''
        nextpept = pepchain+'_'+pepnum2+'_'+pepnum3+'_PEPT'
        if nextpept not in peptide_frags.keys():
            try: 
                seg[res2+'_CTC'] = residues[res2]
                entry = [x[1] for x in residues[res2]]
                name = res2[0]+'_'+aa_dict[residues[res2][0][3]]+res2[1:]
                entry.insert(0,name+'_CTC')
                fA.append(entry) 
                del residues[res2]
                segs.append(seg)
                seg = {}
            # assume terminal peptide
            except:
                continue 
   
    for free in sorted(residues.keys()):
        entry = [x[1] for x in residues[free]]
        name = free[0]+'_'+aa_dict[residues[free][0][3]]+free[1:]
        entry.insert(0,name+'_FREE')
        fA.append(entry)
        del residues[free]
     

    fB = []
    #assumes one ligand
    entry = [x[1] for x in ligand]
    name = 'LIG_'+ligand[0][3]+'_'+ligand[0][5]
    entry.insert(0,name)
    fB.append(entry)

    #assumes solvent goes in c
    fC = []
    for solv in sorted(solvent.keys()):
        entry = [x[1] for x in solvent[solv]]
        name = 'SOLV_'+solvent[solv][0][3]+'_'+solvent[solv][0][5]
        entry.insert(0,name)
        fC.append(entry) 


    cmd.show('sticks') 
    color_frags(fA) 
    color_frags(fB, [0.5,0.5,0.5]) 
    color_frags(fC, [0.0,0.0,1.0]) 
    try: os.mkdir('fsapt/')
    except: print 'Overwriting old fA.dat and fB.dat'
    write_frag_file(fA,'fA')    
    write_frag_file(fB,'fB')
    write_input(protein_for_inp, ligand_for_inp, protein_charge, ligand_charge, solvent_for_inp, file_name)
    
    

#        for reskey in sorted(residues.keys()):
#            reschain = reskey[0]
#            resnum = reskey[1:]

cmd.extend("chop",chop)

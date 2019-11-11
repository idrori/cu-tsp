# Iddo Drori, Columbia University
import numpy as np
import pickle
import Bio.PDB
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.Polypeptide import PPBuilder
from PPBuilderPlus import PPBuilderPlus
from scipy.spatial import distance_matrix
import matplotlib
from matplotlib import pyplot as plt
matplotlib.pyplot.switch_backend('agg')
import DSSP

drive_prefix = '/Volumes/My Passport for Mac/tsp/'
pdb_dir = 'data/pdb/'
pdb_prefix = 'pdb'
pdb_suffix = '.ent'
clean_pdb_suffix = '.clean.pdb'
dssp_dir = 'data/dssp/'
dssp_suffix = '.dssp'
png_dir = 'png/'


from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from pdb import set_trace

embedding1 = {
'A': [1.8, 0.67, 6, 0, 0, 4], 
'C': [2.5, 0.38, 7, 0, 0, 9],
'D': [-3.5, -1.2, 9, -1, 0, 6],
'E': [-3.5, -0.75, 10, -1, 0, 5],
'F': [2.8, -0.57, 12, 0, 1, 6],
'G': [-0.4, 0.48, 5, 0, 0, 6],
'H': [-3.2, 0.64, 11, 1, 1, 8],
'I': [4.5, 1.9, 9, 0, 0, 4],
'K': [-3.9, -0.57, 10, 1, 0, 5],
'L': [3.8, 1.9, 9, 0, 0, 4],
'M': [1.9, 2.4, 9, 0, 0, 5],
'N': [-3.5, -0.6, 9, 0, 0, 6],
'P': [-1.6, 1.2, 8, 0, 0, 7],
'Q': [-3.5, -0.22, 10, 0, 0, 5],
'R': [-4.5, -2.1, 12, 1, 0, 5],
'S': [0.8, 0.01, 7, 0, 0, 4],
'T': [-0.7, 0.52, 8, 0, 0, 5],
'V': [4.2, 1.5, 8, 0, 0, 4],
'W': [-0.9, 2.6, 15, 0, 1, 11],
'Y': [-0.13, 1.6, 13, 0, 1, 7],
# https://www.bioinformatics.org/sms2/iupac.html
'B': [-3.5, -0.9, 9, 0.5, 0, 6], # mean of D and N
'Z': [-3.5, -0.485, 10, -0.5, 0, 5], # mean of Q and E
'X': [4.5, 2.6, 15, 1, 1, 11] # max
}

embedding2 = {
'A': [0.748392, 0.1432352, -1.388808078, -0.09796, -0.487339717, -0.97169], 
'C': [0.983156, -0.0868934, -0.974238503, -0.09796, -0.487339717, 1.72744], 
'D': [-1.0291, -1.3406976, -0.145099351, -2.05714, -0.487339717, 0.107965],
'E': [-1.0291, -0.9836014, 0.269470224, -2.05714, -0.487339717, -0.43186], 
'F': [1.083768, -0.840763, 1.098609375, -0.09796, 1.949358869, 0.107965],
'G': [0.010564, -0.0075387, -1.803377654, -0.09796, -0.487339717, 0.107965], 
'H': [-0.92849, 0.1194288, 0.6840398, 1.86122, 1.949358869, 1.187615],
'I': [1.653908, 1.119298, -0.145099351, -0.09796, -0.487339717, -0.97169],
'K': [-1.16325, -0.840763, 0.269470224, 1.86122, -0.487339717, -0.43186],
'L': [1.419145, 1.119298, -0.145099351, -0.09796, -0.487339717, -0.97169],
'M': [0.78193, 1.5160715, -0.145099351, -0.09796, -0.487339717, -0.43186],
'N': [-1.0291, -0.8645694, -0.145099351, -0.09796, -0.487339717, 0.107965], 
'P': [-0.39189, 0.5638151, -0.559668927, -0.09796, -0.487339717, 0.64779],
'Q': [-1.0291, -0.5630216, 0.269470224, -0.09796, -0.487339717, -0.43186],
'R': [-1.36448, -2.0548898, 1.098609375, 1.86122, -0.487339717, -0.43186],
'S': [-0.12359, -0.3805058, -0.974238503, -0.09796, -0.487339717, -0.97169],
'T': [-0.09005, 0.0242032, -0.559668927, -0.09796, -0.487339717, -0.43186],
'V': [1.553295, 0.8018792, -0.559668927, -0.09796, -0.487339717, -0.97169],
'W': [-0.15712, 1.6747808, 2.342318102, -0.09796, 1.949358869, 2.80709],
'Y': [0.101116, 0.8812339, 1.513178951, -0.09796, 1.949358869, 0.64779],
# https://www.bioinformatics.org/sms2/iupac.html
'B': [-1.0291, -1.1026335,  -0.145099351, -1.07755, -0.487339717, 0.107965], # mean of D and N
'Z': [-1.0291, -0.7733115, 0.269470224, -1.07755, -0.487339717, -0.43186], # mean of Q and E
'X': [1.553295, 1.6747808, 2.342318102, 1.86122, 1.949358869, 2.80709] # max
}

one_hot = dict()
alpha_size = len(embedding1.keys()) + 1
for i, key in enumerate(embedding1.keys()):
    one_hot[key] = np.zeros(alpha_size)
    one_hot[key][i] = 1
one_hot['X'] = np.zeros(alpha_size)

N = 10000
def process_msa(fp, start_index_pnet_aa, end_index_pnet_aa):
    homologs = list()
    with open(fp, 'r') as f:
        buf = list()
        i = 0
        for line in f:    
            if i <= N:
                if line.startswith('>'): 
                    i = i + 1
                    if len(buf) != 0:
                        homolog = "".join(buf) 
                        homolog = "".join(x for x in homolog if not x.islower())
                        homolog = homolog[start_index_pnet_aa:end_index_pnet_aa]
                        homologs.append(homolog)
                    buf = list()
                else:
                    buf.append(line.rstrip())
            else:
                break            
#                if line.startswith('>'):
#                    buf = list()
#                else:
#                    buf.append(line.rstrip())
#    aa = "".join(buf[start_index_pnet_aa:end_index_pnet_aa])
    return homologs#, aa

# Equivalent to doing a matrix multiplication of a n x k matrix with a
# k x n matrix, but where each element of the matrix is in the field
# of size embedding dimension.
def get_covariance(homologs, embedding_dict):
    homologs = [list(x) for x in homologs]
    for (i, y) in enumerate(homologs):
        for (j, x) in enumerate(y):
            try:
                homologs[i][j] = embedding_dict[x]
            except:
                homologs[i][j] = embedding_dict['X']
    homologs = np.array(homologs)
    try:
        covariance = np.tensordot(homologs.transpose([1, 0, 2]), homologs, axes=([1, 2], [0, 2]))
    except:
        covariance = 0 # no homologs
    return covariance

# requires downloading all pdb files using PDBList.py
# get all aa chains
def get_pdb_amino_acid_sequences(pdb_path):
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path[:-4], pdb_path)
    ppb = PPBuilder()
    pdb_aas = []
    for pp in ppb.build_peptides(structure): 
        pdb_aa = str(pp.get_sequence())
        pdb_aas.append(pdb_aa)
    return pdb_aas

def get_calpha_distance_matrix(pdb_path, chain_index):
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path[:-4], pdb_path)
    A = []
    ppb = PPBuilder()
    pdb_aas = []
    model = ppb.build_peptides(structure)
    chain = model[chain_index]
    try:
    	ca_list = chain.get_ca_list()
    except:
    	print('get_calpha_distance_matrix chain.get_ca_list() exception', chain_index)
    	return A
    for ca in ca_list:
        coord = ca.get_coord()
        A.append(np.asarray(coord))
    D = distance_matrix(A,A)
    print('calpha_distance', D[10:15,10:15])
    return D

def get_cbeta_distance_matrix(pdb_path, chain_index):
    print("BEFORE")
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path[:-4], pdb_path)
    A = []
    nobetas_glycine = []
    ppb = PPBuilderPlus()
    pdb_aas = []
    model = ppb.build_peptides(structure)
    chain = model[chain_index]
    try:
    	cb_list = chain.get_cb_list()
    except:
    	print('get_cbeta_distance_matrix chain.get_cb_list() exception', chain_index)
    	return A
    i = 0
    for cb in cb_list:
    	if cb == None:
            nobetas_glycine.append(i)
            coord = np.asarray([0.0, 0.0, 0.0])
    	elif cb != None:
            coord = cb.get_coord()
    	A.append(np.asarray(coord))
    	i = i + 1
    D = distance_matrix(A,A)
    M = np.ones(np.shape(D))
    for i in range(len(nobetas_glycine)):
    	D[:,i] = None
    	D[i,:] = None
    	M[:,i] = 0
    	M[i,:] = 0	
    #print('cbeta_distance', D[10:15,10:15])
    print("AFTER")
    return D, M
    
def get_coords(pdb_path, chain_index):
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path[:-4], pdb_path)
    A = []
    ppb = PPBuilder()
    pdb_aas = []
    model = ppb.build_peptides(structure)
    chain = model[chain_index]
    try:
    	ca_list = chain.get_ca_list()
    except:
    	print('get_coords chain.get_ca_list() exception', chain_index)
    	return A
    for ca in ca_list:
        coord = ca.get_coord()
        A.append(np.asarray(coord))
    return A

def get_pdb_torsion_angles(pdb_path, chain_index):
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path[:-4], pdb_path)
    A = []
    ppb = PPBuilder()
    pdb_aas = []
    model = ppb.build_peptides(structure)
    chain = model[chain_index]
    phi_psi_list = chain.get_phi_psi_list()
    return [x[0] for x in phi_psi_list], [x[1] for x in phi_psi_list]

# requires downloading all dssp files and using DSSP.py
def get_dssp_amino_acid_sequences(pdb_path, dssp_path):
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path[:-4], pdb_path)
    model = structure[0]
    dssp = DSSP.DSSP(model, dssp_path)
    aa = ''
    aas = []
    prev_chain_id = ''
    for key in list(dssp.keys()):
        chain_id = key[0]
        if chain_id != prev_chain_id:
            prev_chain_id = chain_id
            if aa != '':
                aas.append(aa)
            aa = ''
        #if chain_id == letters[chain_index]:#'A': # first chain
        aa += dssp[key][1]
    return aas

# concat all aa chains and their corresponding secondary structures
def get_dssp_amino_acid_sequence(pdb_path, dssp_path):
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path[:-4], pdb_path)
    model = structure[0]
    dssp = DSSP.DSSP(model, dssp_path)
    aa = ''
    for key in list(dssp.keys()):
        #chain_id = key[0]
        #if chain_id == letters[chain_index]:#'A': # first chain
        aa += dssp[key][1]
    return aa

# concat all aa chains and their corresponding secondary structures
def get_dssp_secondary_structure(pdb_path, dssp_path):
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path[:-4], pdb_path)
    model = structure[0]
    dssp = DSSP.DSSP(model, dssp_path)
    q8 = ''
    for key in list(dssp.keys()):
        #chain_id = key[0]
        #if chain_id == letters[chain_index]:#'A': # first chain
        #aa += dssp[key][1]
        q8 += dssp[key][2]
    return q8

def get_dssp_absolute_surface_area(pdb_path, dssp_path):
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path[:-4], pdb_path)
    model = structure[0]
    dssp = DSSP.DSSP(model, dssp_path)
    asa = ''
    #letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' # chain_id is in letters
    for key in list(dssp.keys()):
        #chain_id = key[0]
        #if chain_id == letters[chain_index]:#'A': # first chain
        asa += dssp[key][3]
    return asa

def get_dssp_torsion_angles(pdb_path, dssp_path):
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path[:-4], pdb_path)
    model = structure[0]
    dssp = DSSP.DSSP(model, dssp_path)
    phi = ''
    psi = ''
    for key in list(dssp.keys()):
        phi += dssp[key][4]
        psi += dssp[key][5]
    return phi, psi
    
def get_pdb_bio_amino_acid_sequences(pdb_path):
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path[:-4], pdb_path)
    pdb_aas = []
    for model in structure:
    	for chain in model:
    		pdb_aa = ''
    		for residue in chain:
        		if Bio.PDB.Polypeptide.is_aa(residue):
        			residue_name = residue.get_resname()
        			try:
        				residue_letter = Bio.PDB.Polypeptide.three_to_one(residue_name)
        			except:
        				print('residue_letter = X')
        				residue_letter = 'X'
        			pdb_aa += residue_letter	
    		pdb_aas.append(pdb_aa)
    return pdb_aas

def get_calpha_distance_matrix2(pdb_path, chain_index):
    print(chain_index)
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path[:-4], pdb_path)
    A = []
    for model in structure:
    	print('len(model)', len(model))
    	for index, chain in enumerate(model):
        	if index == chain_index:
        		for residue in chain:
        			try:
        				coord = residue['CA'].get_coord()
        				A.append(np.asarray(coord))
        			except:
        				continue
        		D = distance_matrix(A,A)
        		return D

def get_coords2(pdb_path, chain_index):
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path[:-4], pdb_path)
    A = []
    for model in structure:
    	for index, chain in enumerate(model):
        	if index == chain_index:
        		for residue in chain:
        			try:
        				coord = residue['CA'].get_coord()
        				A.append(np.asarray(coord))
        			except:
        				continue
        		return A

# since there are proteins with many homologs
# do not compute multiplication A^TA directly: expensive, condition number is squared
# use G = A^TA is symmetric pos def, so G = V(S^TS)V^T
# therefore compute G efficiently


SKIPLINES = 33
#SKIPLINES = 48 # 48 for casp13, 33 otherwise.
#filename = 'casp' + '12' + '/' + 'testing'#'training_100'
#infile = open(filename, 'r')
#lines = infile.readlines()
#infile.close()

"""
outfile = open('training_100_1', 'w')
for (i, line) in enumerate(lines):
	if i % (SKIPLINES * 1000) == 0:
		fold = str(int(i / (SKIPLINES * 1000) + 1))
		print("Fold {}".format(fold))
		outfile.close()
		outfile = open('training_100_{}'.format(fold), 'w')
	outfile.write(line)
exit()
"""
#start = 0
#stop = int(len(lines)/SKIPLINES)

NAA = 21
index = 1

import os.path
# using 5 x 4TB fast external drives for uncompressed raw data
driveprefix = '/Volumes/My Passport for Mac'
msapath = 'msa'
suffix = '.a2m'
def doAllMSA(pdbline, start_index_pnet_aa, end_index_pnet_aa):
	prefix = pdbline[0:2].upper()
	print('prefix',prefix)
	# map first letter of pdb to external drive, which get their names based on connection order
	pdb_first_letter = prefix[0:1]
	if pdb_first_letter == '1':
		drivesuffix = ' 4'
	elif pdb_first_letter == '2':
		drivesuffix = ' 3'
	elif pdb_first_letter == '3':
		drivesuffix = ''
	elif pdb_first_letter == '4':
		drivesuffix = ' 1'
	elif pdb_first_letter == '5' or pdb_first_letter == '6' or pdb_first_letter == '7' or pdb_first_letter == '8' or pdb_first_letter == '9':
		drivesuffix = ' 2'
	drivepath = driveprefix + drivesuffix
	fname = drivepath + '/' + msapath + '/' + prefix + '/' + pdbline.upper() + suffix
	print('fname',fname)
	print('os.path.isfile(fname)',os.path.isfile(fname))
	print('start_index_pnet_aa',start_index_pnet_aa)
	print('end_index_pnet_aa',end_index_pnet_aa)
	try:
		homologs = process_msa(fname, start_index_pnet_aa, end_index_pnet_aa)
		print('len(homologs)', len(homologs))
		#print('aa', aa)
		msacov_onehot = get_covariance(homologs, one_hot)
		msacov_embedding = get_covariance(homologs, embedding2)
		if (len(homologs) == 0): # no homologs in data
			msacov_onehot = np.zeros((n,n))
			msacov_embedding = np.zeros((n,n))			
		msa_flag = True
	except:
		n = end_index_pnet_aa - start_index_pnet_aa
		msacov_onehot = np.zeros((n,n))
		msacov_embedding = np.zeros((n,n))
		msa_flag = False
	print('msa_flag', msa_flag)
	return msacov_onehot, msacov_embedding, msa_flag

# output: pickled python dictionary of {"protein_id": {"length": 1, "aa": ni, "ss": ni, "dcalpha": ni*ni, "x": ni, "y": ni, "z": ni, "psi": ni, "phi": ni, "pssm": 21xni}}

def rama(pdb, phi, psi): # input in radians
	plt.scatter(phi, psi, marker=".")
	plt.xlim(-np.pi, np.pi)
	plt.xlabel("phi")
	plt.ylabel("psi")
	plt.ylim(-np.pi, np.pi)
	plt.savefig(png_dir + pdb + '.rama.png', bbox_inches='tight')

def get_dihedral(coords1, coords2, coords3, coords4):
    a1 = coords2 - coords1
    a2 = coords3 - coords2
    a3 = coords4 - coords3
    v1 = np.cross(a1, a2)
    v1 = v1 / (v1 * v1).sum(-1)**0.5
    v2 = np.cross(a2, a3)
    v2 = v2 / (v2 * v2).sum(-1)**0.5
    porm = np.sign((v1 * a3).sum(-1))
    rad = np.arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
    if not porm == 0:
        rad = rad * porm
    return rad

def separate_coords(coords, pos):
    res = []
    for i in range(len(coords[0])):
        if i % 3 == pos:
            res.append([coords[j][i] for j in range(3)])
    return np.array(res)
    
def get_pnet_angles(coords):
    coords_nterm = separate_coords(coords, 0)
    coords_calpha = separate_coords(coords, 1)
    coords_cterm = separate_coords(coords, 2)
    ph_angle_dists, ps_angle_dists = [], []
    phi, psi = [0.0], []
    for i in range(len(coords_calpha)):
        if i > 0:
            phi.append(get_dihedral(coords_cterm[i-1], coords_nterm[i], coords_calpha[i], coords_cterm[i])) # my_calc
        if i < len(coords_calpha) - 1: 
            psi.append(get_dihedral(coords_nterm[i], coords_calpha[i], coords_cterm[i], coords_nterm[i+1])) # my_calc
    psi.append(0)
    return phi, psi

for fold in range(1, 105): # Change 2 back to 105 to training
	data = {}
	#indices = []
	#pdbs = []
	#length_aas = []
	#pdb_aas = []
	#dcalphas = []
	#coords = []
	#q8s = []
	#psis = []
	#phis = []
	#pssms = []
	k = 1000 #for training, 223 for validation, 30 for testing
	start = 0
	stop = k
	print("Fold: {}".format(fold))
	lines = open(drive_prefix+'training_100_{}'.format(fold)).readlines()
	#lines = open('validation').readlines()
	lines = open('testing_casp13.proteinnet').readlines()
	#outfile = open('training_100msa_'+str(fold)+'.csv', 'w')
	#outfile = open('validation.csv', 'w')
	for i in range(start, stop):
		#if i == start+2:
		#	exit()
		d = SKIPLINES*i
		s = 0 # 0 for training, 3 for validation, 5 for testing
		pdbline = lines[1+d].lower().strip()[0+s:] # 5+s: for testing only
		pdb = lines[1+d].lower().strip()[0+s:4+s]
		print(i, pdb)
		pdb_path = drive_prefix + pdb_dir + pdb_prefix + pdb + pdb_suffix
		#pdb_path = drive_prefix + pdb_dir + pdb_prefix + pdb + pdb_suffix
		dssp_path = drive_prefix + dssp_dir + pdb + dssp_suffix
		#dssp_path = '"' + drive_prefix + dssp_dir + pdb + dssp_suffix + '"'
		clean_pdb_path = '"' + drive_prefix + pdb_dir + pdb_prefix + pdb + clean_pdb_suffix + '"'
	
		# protein data bank has pdb's with missing pdb format files,
		# as well as pdb's which have been superseded with other id's
		print(pdb_path)
		print(dssp_path)
		try:
			pdb_aas = get_pdb_amino_acid_sequences(pdb_path)
		except:
			print('pdb format/id exception', pdb)
			pdb_aas = []
		try:
			dssp_aas = get_dssp_amino_acid_sequences(pdb_path, dssp_path)
		except:
			print('dssp exception', pdb)
			dssp_aas = []

		pnet_aa = str(lines[3+d].strip())

		print('pdb_aas', pdb_aas, len(pdb_aas))
		print('dssp_aas', dssp_aas, len(dssp_aas))		
		print('pnet_aa', pnet_aa)
		#exit()

		foundmatch = False
		sorted_pdb_aas = sorted(enumerate(pdb_aas), key=lambda x: -len(x[1]))
		#for chain_index in range(len(pdb_aas)): # for each aa chain in pdb
		for chain_index, _ in sorted_pdb_aas:
			if foundmatch == False:
				pdb_aa = pdb_aas[chain_index]
				print('pdb_aa', chain_index, pdb_aa)
								
				if pdb_aa and pdb_aa in pnet_aa: # not empty and pdb aa sequence substring of pnet aa sequence
					print('pdb_aa in pnet_aa')
					# amino acid (aa) sequence
					length_aa = len(pdb_aa)
			
					# pdb aa is a substring of pnet aa which allows using corresponding subset of pssm features
					start_index_pnet_aa = pnet_aa.find(pdb_aa)
					end_index_pnet_aa = start_index_pnet_aa + length_aa
					sub_pnet_aa = pnet_aa[start_index_pnet_aa:end_index_pnet_aa]

					# pssm features: 21 x sequence length
					pnet_pssm = lines[5+d:5+NAA+d]
					sub_pnet_pssms = []
					for p in range(NAA): # for each amino acid type
						pnet_pssm_p = pnet_pssm[p].split()
						sub_pnet_pssm_p = pnet_pssm_p[start_index_pnet_aa:end_index_pnet_aa] # length of aa sequence
						sub_pnet_pssms.append(sub_pnet_pssm_p)
			
					pnet_coord = lines[5+NAA+d+1:5+NAA+d+4]
					for i in range(0,3):
						pnet_coord[i] = pnet_coord[i].split("\t")
						for j in range(0,len(pnet_coord[i])):
							pnet_coord[i][j] = float(pnet_coord[i][j])
					#pnet_phi, pnet_psi = get_pnet_angles(pnet_coord)
					#sub_pnet_phi = pnet_phi[start_index_pnet_aa:end_index_pnet_aa]
					#sub_pnet_psi = pnet_psi[start_index_pnet_aa:end_index_pnet_aa]
					
					# distance matrices
					dcalpha = get_calpha_distance_matrix(pdb_path, chain_index)
					dcbeta, dcbeta_mask = get_cbeta_distance_matrix(pdb_path, chain_index)
					#dcbeta = None

					# 3d coordinates
					coord = get_coords(pdb_path, chain_index)		
								
					phi, psi = get_pdb_torsion_angles(pdb_path, chain_index) # in radians
					in_units = 1#180/np.pi					
					for i in range(0, len(phi)): # convert to degrees, first is None
						if phi[i] != None:
							phi[i] = phi[i] * in_units
						else:
							phi[i] = 0
					for i in range(0, len(psi)): # convert to degrees, last is None
						if psi[i] != None:
							psi[i] = psi[i] * in_units
						else:
							psi[i] = 0				
					#rama(pdb, phi, psi)
					#print(phi[10]-sub_pnet_phi[10]) # difference is 10^-8 make sure this is not significant
					#exit()
	
					try:
						dssp_aa = get_dssp_amino_acid_sequence(pdb_path, dssp_path)
						dssp_q8 = get_dssp_secondary_structure(pdb_path, dssp_path)
						#dssp_asa = get_dssp_absolute_surface_area(pdb_path, dssp_path)
						#dssp_phi, dssp_psi = get_dssp_torsion_angles(pdb_path, dssp_path)
					except:
						print('dssp exception', dssp_path)
						dssp_aa = ''
						dssp_q8 = ''
						#dssp_asa = ''
						#dssp_phi_psi = ''
					q8 = ''
					asa = ''
					if (pdb_aa in dssp_aa):
						# pdb aa is a substring of dssp_aa which allows using corresponding subset of secondary structures
						start_index_dssp_aa = dssp_aa.find(pdb_aa)
						end_index_dssp_aa = start_index_dssp_aa + length_aa
						sub_dssp_aa = dssp_aa[start_index_dssp_aa:end_index_dssp_aa]
						q8 = dssp_q8[start_index_dssp_aa:end_index_dssp_aa]
						#asa = dssp_asa[start_index_dssp_aa:end_index_dssp_aa]
						#dssp_phi = dssp_phi[start_index_dssp_aa:end_index_dssp_aa]
						#dssp_psi = dssp_psi[start_index_dssp_aa:end_index_dssp_aa]
						
					# process msa features of protein with corresponding chain
					# to generate 2 n x n co-variance matrices: for 1-hot encoding and embedding
					#print('pdbline',pdbline)
					#msacov_onehot, msacov_embedding, msa_flag = doAllMSA(pdbline, start_index_pnet_aa, end_index_pnet_aa)
					# find file in correct drive based on protein prefix
					# take start and end indices
					# compute covariance matrix
					#print('pdb',pdb)
					#print('chain_index', chain_index)
					#print('pdb_aa',pdb_aa)
					#fp = '' # concat protein name and chain
					#homologs, pdb = process_msa(fp)
					#msacov_onehot = get_covariance(homologs, one_hot)
					#cov_embedding1 = get_covariance(homologs, embedding1)
					#msacov_embedding = get_covariance(homologs, embedding2)					
						
					# verify equal lengths
					#print('aa len',length_aa,len(pnet_aa))
					#print('dcalpha len',length_aa,dcalpha.shape[0],dcalpha.shape[1])
					#print('phi len',length_aa,len(phi))
					#print('psi len',length_aa,len(psi))
					#print('len(dcalpha)', len(dcalpha))
					if (len(dcalpha) > 0) and (length_aa <= len(pnet_aa)) and (length_aa == dcalpha.shape[0]) and (length_aa == dcalpha.shape[1]) and (length_aa == len(phi)) and (length_aa == len(psi)) and (length_aa == len(q8)): 
						msacov_onehot, msacov_embedding, msa_flag = doAllMSA(pdbline, start_index_pnet_aa, end_index_pnet_aa)
						#if msa_flag:
						if True:
							#if (q8 == '' or length_aa != len(q8)):
							#	print('pdb',pdb,'chain_index',chain_index,'length_aa',length_aa,'len(dssp_aa)',len(dssp_aa),'len(q8)',len(q8),'q8',q8)
							#	print('dssp_ss',dssp_ss)
							#	print('ss     ',ss)
							#	print('dssp_aa',dssp_aa)
							#	print('pdb_aa ',pdb_aa)
							#	exit()
							foundmatch = True
							print('match found')
							record = {}										
							record['length'] = length_aa 					# scalar
							record['chain'] = chain_index 					# scalar
							# data record for each protein
							# inputs:
							record['aa'] = pdb_aa 							# n x 1
							record['pssm'] = sub_pnet_pssms 				# n x 21					
							record['msacovonehot'] = msacov_onehot 			# n x n
							record['msacovembedding'] = msacov_embedding 	# n x n
							# for both
							record['ss'] = q8 								# n x 1
							# outputs
							record['dcalpha'] = dcalpha 					# n x n
							record['dcbeta'] = dcbeta 					# n x n
							record['dcbeta_mask'] = dcbeta_mask 					# n x n
							record['coords'] = coord 						# n x 3
							record['phi'] = phi 							# n x 1
							record['psi'] = psi 							# n x 1
							data[pdb] = record
							print(record['length'], len(record['aa']), 'record[aa]', record['aa'])
							#length_aas.append(length_aa)			
							#pdbs.append(pdb)
							#pdb_aas.append(pdb_aa)
							#pssms.append(sub_pnet_pssms)
							#dcalphas.append(dcalpha)
							#coords.append(coord)
							#q8s.append(q8)
							#psis.append(psi)
							#phis.append(phi)
							#indices.append(index)
					
						# png image of distance matrix (for visualization)
						#plt.imshow(dcalpha, cmap='plasma')
						#plt.savefig(png_dir + pdb + '.dcalpha.png', bbox_inches='tight')
						#plt.imshow(msacov_onehot, cmap='plasma')
						#plt.savefig(png_dir + pdb + '.msacov_onehot.png', bbox_inches='tight')
						#plt.imshow(msacov_embedding, cmap='plasma')
						#plt.savefig(png_dir + pdb + '.msacov_embedding.png', bbox_inches='tight')
						# text file (for reading and visualization)
							print(index, pdb, chain_index, length_aa, pdb_aa)
						#outfile.write(str(index)+','+str(pdb)+','+str(chain_index)+','+str(length_aa)+','+str(pdb_aa)+'\n')
						#outfile.flush()
							print(index)
							index = index + 1	

	#outfile.close()
	
	#data = [indices, pdbs, length_aas, pdb_aas, q8s, dcalphas, coords, psis, phis, pssms]
	#datafile = 'cuprotein_'+str(fold)+'.pkl'	
	#datafile = 'cuprotein_validation.pkl'
	datafile = 'cuprotein_testing_casp13.pkl'
	with open(datafile, 'wb') as dataoutfile:
		pickle.dump(data, dataoutfile, pickle.HIGHEST_PROTOCOL)

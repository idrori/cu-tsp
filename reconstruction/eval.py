import numpy as np
import argparse
import math
import torch
import torch.nn.functional as F
import Bio.PDB
import pickle
import subprocess
import PeptideBuilder
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix

from Bio.PDB import *
from Bio.PDB import PDBParser
from Bio.PDB import Superimposer
from Bio.PDB.Atom import *
from Bio.PDB.Residue import *
from Bio.PDB.Chain import *
from Bio.PDB.Model import *
from Bio.PDB.Structure import *
from Bio.PDB.Vector import *
from Bio.PDB.Entity import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.PDBList import PDBList
from pyrosetta import *
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.utility import vector1_core_id_AtomID, vector1_numeric_xyzVector_double_t
from pyrosetta.rosetta import numeric
init()
#from recon_tester import recover_coords, align
from pdb import set_trace

from skimage.restoration import denoise_tv_chambolle, estimate_sigma

import random

from sklearn import manifold
def mds(D):
	n = D.shape[0]
	seed = np.random.RandomState(seed=3)
	mds = manifold.MDS(n_components=3, max_iter=3000, eps=1e-9, random_state=seed, dissimilarity="precomputed", n_jobs=1)
	Xtag = mds.fit(D).embedding_
	return Xtag

import cvxpy as cp
def sdp(D):
    lam = 1
    n = D.shape[0]
    G = cp.Variable((n, n), symmetric=True)
    eta = cp.Variable((n, n))
    constraints = [G >> 0] 
    constraints += [eta[i, j] >= 0 for i in range(n) for j in range(n)]
    reg = lam * cp.norm(eta, p=1)
    Gii = cp.diag(G)
    hstacked = cp.hstack([cp.reshape(Gii, (n, 1)) for i in range(n)])
    vstacked = cp.vstack([cp.reshape(Gii, (1, n)) for i in range(n)])
    obj = reg + 0.5 * cp.sum((hstacked + vstacked - 2 * G + eta - D**2)**2)
    #obj = 0.5 * cp.sum((hstacked + vstacked - 2 * G - D**2))
    prob = cp.Problem(cp.Minimize(obj), constraints)
    prob.solve(verbose=False)
    #print("The optimal value is", prob.value)
    #print("A solution G is: {}".format(G.value))
    u, s, v = np.linalg.svd(G.value)
    u, s = np.real(u), np.real(s)
    dim = 3
    Xtag = np.dot(u, np.diag(np.sqrt(s)))[:,0:dim]
    return Xtag

def sdr(D):
	lam = 1
	n = D.shape[0]
	x = -1./(n + np.sqrt(n))
	y = -1./np.sqrt(n)
	V = np.vstack([y * np.ones((1, n-1)), x * np.ones((n-1,n-1)) + np.identity(n-1)])
	e = np.ones((n, 1))
	eT = np.transpose(e)
	G = cp.Variable((n-1, n-1), symmetric=True)
	B = V * G * np.transpose(V)
	diagB = cp.atoms.affine.diag.diag(B)
	diagB = cp.atoms.affine.reshape.reshape(diagB, (n,1))
	diagBT = diagB.T
	E = diagB * eT + e * diagBT - 2 * B
	constraints = [G >> 0]
	obj = cp.trace(G) - lam * cp.atoms.norm(E - D, "fro")
	prob = cp.Problem(cp.Maximize(obj), constraints)
	prob.solve(verbose=True)
	U, S, V = np.linalg.svd(B.value)
	dim = 3
	Xtag = (np.sqrt(S)*np.transpose(V))[:,0:dim]
	return Xtag

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

def separate_coords(X_fullatom, pos):
    res = []
    n = X_fullatom.shape[0]
    for i in range(n):
        if i % 3 == pos:
            res.append(X_fullatom[i])
    return np.array(res)

def compute_torsion_angles(X_fullatom):
    X_nterm = separate_coords(X_fullatom, 0)
    X_calpha = separate_coords(X_fullatom, 1)
    X_cterm = separate_coords(X_fullatom, 2)
    phis = [0]
    psis = []
    n = len(X_calpha)
    for i in range(n):
        if i > 0:
            phis.append(get_dihedral(X_cterm[i-1], X_nterm[i], X_calpha[i], X_cterm[i]))
        if i < n - 1: 
            psis.append(get_dihedral(X_nterm[i], X_calpha[i], X_cterm[i], X_nterm[i+1]))
    psis.append(0)
    return phis, psis

def load_data(path, is_pred=False, with_Distance=True, with_Torsion=False):
    print(path)
    data = pickle.load(open(path, 'rb'))
    if path == 'cgan_100k.pkl':
        for i in range(len(data)):
            print(data[i], data[i].shape)
            d = np.reshape(data[i], (128, 128))
            print(d.shape)
            saveD(str(i), d, str(i))
        exit()
    pdbs = sorted(list(data.keys()))
    out= dict()
    if not is_pred:
        out['lengths'] = []
        out['chains'] = []
        out['aas'] = []
        out['q8s'] = []
        out['pssms'] = []

    out['pdbs'] = pdbs
    out['dcalphas'] = []
    out['coords'] = []
    out['psis'] = []
    out['phis'] = []
    for p in pdbs:
        if not is_pred:
            out['lengths'].append(data[p]['length'])
            out['chains'].append(data[p]['chain'])
            out['aas'].append(data[p]['aa'])
            out['q8s'].append(data[p]['ss'])
            out['pssms'].append(data[p]['pssm'])
            out['coords'].append(data[p]['coords'])
        if 'model100k_' in path:
            if with_Distance:
                out['dcalphas'].append(data[p]['dcalpha'])
            if with_Torsion:
                out['psis'].append(data[p]['psi'])
                out['phis'].append(data[p]['phi'])
        else:
            if with_Distance:
                out['dcalphas'].append(data[p]['dcalpha'])
            if with_Torsion:
                out['psis'].append(data[p]['psi'])
                out['phis'].append(data[p]['phi'])
    return out

def calc_pairwise_distances(chain_a, chain_b, use_gpu):
    distance_matrix = torch.Tensor(chain_a.size()[0], chain_b.size()[0]).type(torch.float)
    # add small epsilon to avoid boundary issues
    epsilon = 10 ** (-4) * torch.ones(chain_a.size(0), chain_b.size(0))
    if use_gpu:
        distance_matrix = distance_matrix.cuda()
        epsilon = epsilon.cuda()

    for i, row in enumerate(chain_a.split(1)):
        distance_matrix[i] = torch.sum((row.expand_as(chain_b) - chain_b) ** 2, 1).view(1, -1)

    return torch.sqrt(distance_matrix + epsilon)

def calc_drmsd(chain_a, chain_b, use_gpu=False):
    assert len(chain_a) == len(chain_b)
    chain_a = torch.from_numpy(chain_a)
    chain_b = torch.from_numpy(chain_b)
    distance_matrix_a = calc_pairwise_distances(chain_a, chain_a, use_gpu)
    distance_matrix_b = calc_pairwise_distances(chain_b, chain_b, use_gpu)
    return torch.norm(distance_matrix_a - distance_matrix_b, 2) \
            / math.sqrt((len(chain_a) * (len(chain_a) - 1)))

def calc_gdt_ts(a, b):
    N = a.shape[0]
    count1, count2, count4, count8 = 0, 0, 0, 0
    for i in range(N):
        d = np.sum((a[i, :] - b[i, :])**2)**0.5
        if d <= 1:
            count1 += 1
        if d <= 2:
            count2 += 1 
        if d <= 4:
            count4 += 1
        if d <= 8:
            count8 += 1
    gdt_ts = ((count1 + count2 + count4 + count8) / (4*N)) * 100
    return gdt_ts

def set_pdb_coords(pdb, coords, out_path):
    pose = pose_from_pdb(pdb)
    n = coords.shape[0]
    for i in range(1, n + 1):
        x, y, z = coords[i - 1, :]
        out = numeric.xyzVector_double_t(x, y, z)
        pose.residues[i].atom("CA").xyz(out)
    pose.dump_pdb(out_path)

def to_pdb(aa, phi, psi, path):
    phi, psi = phi[1:], psi[:-1]
    struct = PeptideBuilder.make_structure(aa, phi, psi)
    out = Bio.PDB.PDBIO()
    out.set_structure(struct)
    out.save(path)

def print_metrics(p, gt_path, pred_path, a_data, b_data, method):
    gt_pdb_path = '{}_{}.pdb'.format(gt_path[:-4], p)
    pred_pdb_path = '{}_{}_{}.pdb'.format(pred_path[:-4], p, method)
    gt_aa, gt_phi, gt_psi, gt_coords = a_data
    if method == 'torsion':
        b_aa, b_phi, b_psi = b_data
        gt_path = gt_path[:-4]
        pred_path = pred_path[:-4]
        to_pdb(gt_aa, gt_phi, gt_psi, gt_pdb_path)
        to_pdb(b_aa, b_phi, b_psi, pred_pdb_path)
    elif method == 'dcalpha':
        pred_dcalpha = b_data
        pred_coords = recover_coords(pred_dcalpha, method='SDP') 
        #aligned_coords = align(pred_coords, gt_coords)
        aligned_coords = pred_coords

        to_pdb(gt_aa, gt_phi, gt_psi, gt_pdb_path)
        set_pdb_coords(gt_pdb_path, np.array(gt_coords), pred_pdb_path)
    print("-----------{}-{}------------".format(p, method))
    proc = subprocess.Popen('java -jar TMscore.jar {} {}'.format(gt_pdb_path, pred_pdb_path),
                stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.decode('utf-8').split('\n')
    for i in range(10, 18):
        print(out[i].rstrip())
    print("---------------------------".format(p))

def torsion_only(aa, phi, psi):
    phi, psi = phi[1:], psi[:-1]
    recon_struct = PeptideBuilder.make_structure(aa, phi, psi)
    atoms = recon_struct.get_atoms()
    coords = list()
    for atom in atoms:
        coords.append(atom.get_coord())
    return np.array(coords)

d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

d1to3 = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
     'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
     'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
     'A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def saveD(p, D, text):
    """
    plt.cla()
    plt.clf()
    plt.close()   
    plt.imshow(D, cmap='plasma')
    plt.savefig('{}_D_'.format(p)+text+'.png')
    """

def saveScatter(p, X, text):
    n = X.shape[0]
    xs = []
    ys = []
    zs = []
    for j in range(n):
        x,y,z = X[j][0], X[j][1], X[j][2]
        xs.append(x)
        ys.append(y)
        zs.append(z)
    """
    plt.cla()
    plt.clf()
    plt.close()    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.cla()
    ax.scatter(xs,ys,zs, color='blue', label=text)
    ax.set_title(p)
    ax.legend(loc='upper left', fontsize='x-large')
    plt.savefig('{}_X_'.format(p)+text+'.png')
    """

import PeptideBuilder
from PeptideBuilder import Geometry
import Bio.PDB

def viz(p, aa, D, X, text, reflect = False):
    xs = []
    ys = []
    zs = []
    file = open('{}_'.format(p)+text+'.xyz', 'w')
    n = len(aa)
    file.write(str(n) + '\n')
    file.write(text + '_' + p + '\n')
    for j in range(n):
        x,y,z = X[j][0], X[j][1], X[j][2]
        if reflect == True:
            x = -x
        file.write(aa[j] + '\t' + str(x) + '\t' + str(y) + '\t' + str(z) + '\n') 
        xs.append(x)
        ys.append(y)
        zs.append(z)
    file.close()
    saveD(p, D, text)

def reflect(X, axis):
    n = X.shape[0]
    Xtag = np.zeros((n,3))
    for i in range(n):
        Xtag[i][axis] = -X[i][axis]
    return Xtag

import PeptideBuilderCoords

def build_fullatom_from_distance(p, aa, X, text):
    geos = []
    for a in aa:
        geos.append(Geometry.geometry(a))
    structure = PeptideBuilderCoords.make_structure_from_geos_coords(geos, X) 
    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save('{}_'.format(p)+text+'.pdb')

def build_fullatom_from_torsion(p, aa, X, text):
    #phis, psis = compute_torsion_angles(X_fullatom)
    #print(n, len(phis), len(psis), phis)
    #structure = PeptideBuilder.make_structure(aa, phis, psis)
    geos = []
    for a in aa:
        geos.append(Geometry.geometry(a))
    structure = PeptideBuilderCoords.make_structure_from_geos_coords(geos, X) 
    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save('{}_'.format(p)+text+'.pdb')


from Bio.SVDSuperimposer import SVDSuperimposer

def drmsd(X, Xtag):
    n = X.shape[0]
    return (1.0/n * np.sum((X - Xtag)**2))**0.5

def superimpose(X, Xtag):
    try:
        sup = SVDSuperimposer()
        sup.set(X, Xtag)
        sup.run()
        Xtag_super = sup.get_transformed()
        return Xtag_super
    except:
        return np.zeros((X.shape[0],X.shape[1]))

def relative_error(Atag, A):
    try:
        return np.linalg.norm(Atag - A) / np.linalg.norm(A)
    except:
        print('relative error except')
        return 1

def torsion_relative_error(Atag, A, Btag, B):
    try:
        e = 0
        for i in range(len(A)):
            if Atag[i] != None and A[i] != None:
                e = e + abs(Atag[i] - A[i])
        for i in range(len(B)):
            if Btag[i] != None and B[i] != None:
                e = e + abs(Btag[i] - B[i])
        e = e / (len(A) + len(B))
        e = e*np.pi/180
        return e
    except:
        print('torsion relative error except')
        return 1

def triangle_inequality(d1, d2, d3):
    if d1 > d2 and d1 > d3 and d1 > (d2 + d3):
        print(d1, (d2 + d3))
        return False
    elif d2 > d1 and d2 > d3 and d2 > (d1 + d3):
        print(d2, (d1 + d3))
        return False
    elif d3 > d2 and d3 > d1 and d3 > (d1 + d2):
        print(d3, (d1 + d2))
        return False 
    return True

# add probs to pkl

def jitter(D):
    n = D.shape[0]
    mu, sigma = 0, 1
    s = np.random.normal(mu, sigma, n)

def improve_energy(pose):    
    scorefxn = get_fa_scorefxn()
    scorefxn(pose)

def distance_only(distance_pred_path, gt_path):
    distance_pred_data = load_data(distance_pred_path, is_pred=True)
    gt_data = load_data(gt_path, is_pred=False)
    for (i, p) in enumerate(gt_data['pdbs']):
        aa = gt_data['aas'][i]
        print(p, aa)
        n = len(aa)
        D = gt_data['dcalphas'][i]                 

        X = np.reshape(gt_data['coords'][i], (n,3))
        build_fullatom_from_distance(p, aa, X, distance_pred_path[:-4]+'_'+'ground_truth')
        viz(p, aa, D, X, distance_pred_path[:-4]+'_'+'ground_truth')
        print('||distance_matrix(X,X) - D||/||D||', relative_error(distance_matrix(X,X), D))
             
        Xtag = mds(D)
             #Xtagrx = reflect(Xtag,0)
             #Xtagry = reflect(Xtag,1)
             #Xtagrz = reflect(Xtag,2)
             #if relative_error(Xtag, X) > relative_error(Xtagrx,X):
             #    Xtag = Xtagrx
             #elif relative_error(Xtag, X) > relative_error(Xtagry,X):
             #    Xtag = Xtagry
             #elif relative_error(Xtag, X) > relative_error(Xtagrz,X):
             #    Xtag = Xtagrz
        Xtag = superimpose(X, Xtag)
        build_fullatom_from_distance(p, aa, Xtag, distance_pred_path[:-4]+'_'+'ground_truth_reconstructed')
        Dtag = distance_matrix(Xtag, Xtag)
        viz(p, aa, Dtag, Xtag, distance_pred_path[:-4]+'_'+'ground_truth_reconstructed', False)
        print('Xtag=recon(D) Dtag=distance_matrix(Xtag, Xtag) ||Dtag - D||/||D||', relative_error(Dtag, D))
             
        D_pred = distance_pred_data['dcalphas'][i]

#        for k in range(n):
#             a = random.randint(0, n-1) 
#             b = random.randint(0, n-1) 
#             c = random.randint(0, n-1) 
#             if triangle_inequality(D_pred[a][b],D_pred[b][c],D_pred[a][c]) == False:
#                 print(k,a,b,c,D_pred[a][b],D_pred[b][c],D_pred[a][c])
#                 exit()
    
        print('||D_pred - D||/||D||', relative_error(D_pred, D))
        sigma_est = estimate_sigma(D_pred)
        D_pred_denoise = denoise_tv_chambolle(D_pred, weight = sigma_est)
        print('||D_pred_denoise - D||/||D||', relative_error(D_pred_denoise, D))
             
        X_pred_tag = mds(D_pred)
        X_pred_tag = superimpose(X, X_pred_tag)
        build_fullatom_from_distance(p, aa, X_pred_tag, distance_pred_path[:-4]+'_'+'predicted_reconstructed')
        D_pred_tag = distance_matrix(X_pred_tag, X_pred_tag)
        print('||D_pred_tag - D_pred||/||D_pred||', relative_error(D_pred_tag, D_pred))
        print('||D_pred_tag - D||/||D||', relative_error(D_pred_tag, D))
        viz(p, aa, D_pred_tag, X_pred_tag, distance_pred_path[:-4]+'_'+'predicted_reconstructed', False)


def distance_only_comparison(d1, d2, d3, d4, d5, d6, d7, d8, gt):
    distance_pred_data_1 = load_data(d1, is_pred=True)
    distance_pred_data_2 = load_data(d2, is_pred=True)
    distance_pred_data_3 = load_data(d3, is_pred=True)
    distance_pred_data_4 = load_data(d4, is_pred=True)
    distance_pred_data_5 = load_data(d5, is_pred=True)
    distance_pred_data_6 = load_data(d6, is_pred=True)
    distance_pred_data_7 = load_data(d7, is_pred=True)
    distance_pred_data_8 = load_data(d8, is_pred=True)

    gt_data = load_data(gt, is_pred=False)
    
    for (i, p) in enumerate(gt_data['pdbs']):
        aa = gt_data['aas'][i]
        print(p, aa)
        n = len(aa)
        print(p, n)
        D = gt_data['dcalphas'][i] 
        X = np.reshape(gt_data['coords'][i], (n,3))
        build_fullatom_from_distance(p, aa, X, d1[:-4]+'_'+'ground_truth')
        print('||distance_matrix(X,X) - D||/||D||', relative_error(distance_matrix(X,X), D))
             
        Xtag = mds(D)
        Xtag = superimpose(X, Xtag)
        build_fullatom_from_distance(p, aa, Xtag, d1[:-4]+'_'+'ground_truth_reconstructed')
        Dtag = distance_matrix(Xtag, Xtag)
        print('Xtag=recon(D) Dtag=distance_matrix(Xtag, Xtag) ||Dtag - D||/||D||', relative_error(Dtag, D))
             
        D_pred_1 = distance_pred_data_1['dcalphas'][i]
        D_pred_2 = distance_pred_data_2['dcalphas'][i]
        D_pred_3 = distance_pred_data_3['dcalphas'][i]
        D_pred_4 = distance_pred_data_4['dcalphas'][i]
        D_pred_5 = distance_pred_data_5['dcalphas'][i]
        D_pred_6 = distance_pred_data_6['dcalphas'][i]
        D_pred_7 = distance_pred_data_7['dcalphas'][i]
        D_pred_8 = distance_pred_data_8['dcalphas'][i]

        r1 = relative_error(D_pred_1, D)
        r2 = relative_error(D_pred_2, D)
        r3 = relative_error(D_pred_3, D)
        r4 = relative_error(D_pred_4, D)
        r5 = relative_error(D_pred_5, D)
        r6 = relative_error(D_pred_6, D)
        r7 = relative_error(D_pred_7, D)
        r8 = relative_error(D_pred_8, D)

        r = 1.0
        if r1 < r:
             r = r1
             distance_pred_path = d1
             D_pred = D_pred_1
        if r2 < r:
             r = r2
             distance_pred_path = d2
             D_pred = D_pred_2
        if r3 < r:
             r = r3
             distance_pred_path = d3
             D_pred = D_pred_3
        if r4 < r:
             r = r4
             distance_pred_path = d4
             D_pred = D_pred_4
        if r5 < r:
             r = r5
             distance_pred_path = d5
             D_pred = D_pred_5
        if r6 < r:
             r = r6
             distance_pred_path = d6
             D_pred = D_pred_6
        if r7 < r:
             r = r7
             distance_pred_path = d7
             D_pred = D_pred_7
        if r8 < r:
             r = r8
             distance_pred_path = d8
             D_pred = D_pred_8 
        print('BEST ||D_pred - D||/||D||', p, distance_pred_path, relative_error(D_pred, D))
             
        X_pred_tag = mds(D_pred)
                 
        X_pred_tag = superimpose(X, X_pred_tag)
        build_fullatom_from_distance(p, aa, X_pred_tag, distance_pred_path[:-4]+'_'+'predicted_reconstructed')
        D_pred_tag = distance_matrix(X_pred_tag, X_pred_tag)
        print('||D_pred_tag - D_pred||/||D_pred||', relative_error(D_pred_tag, D_pred))
        print('||D_pred_tag - D||/||D||', relative_error(D_pred_tag, D))

def torsion_only_comparison(t1, t2, t3, t4, t5, t6, t7, t8, gt):
    torsion_pred_data_1 = load_data(t1, is_pred=True)
    torsion_pred_data_2 = load_data(t2, is_pred=True)
    torsion_pred_data_3 = load_data(t3, is_pred=True)
    torsion_pred_data_4 = load_data(t4, is_pred=True)
    torsion_pred_data_5 = load_data(t5, is_pred=True)
    torsion_pred_data_6 = load_data(t6, is_pred=True)
    torsion_pred_data_7 = load_data(t7, is_pred=True)
    torsion_pred_data_8 = load_data(t8, is_pred=True)

    gt_data = load_data(gt, is_pred=False)
    
    for (i, p) in enumerate(gt_data['pdbs']):
        S = gt_data['psis'][i] 
        H = gt_data['phis'][i] 

        S_pred_1 = torsion_pred_data_1['psis'][i]
        S_pred_2 = torsion_pred_data_2['psis'][i]
        S_pred_3 = torsion_pred_data_3['psis'][i]
        S_pred_4 = torsion_pred_data_4['psis'][i]
        S_pred_5 = torsion_pred_data_5['psis'][i]
        S_pred_6 = torsion_pred_data_6['psis'][i]
        S_pred_7 = torsion_pred_data_7['psis'][i]
        S_pred_8 = torsion_pred_data_8['psis'][i]

        H_pred_1 = torsion_pred_data_1['phis'][i]
        H_pred_2 = torsion_pred_data_2['phis'][i]
        H_pred_3 = torsion_pred_data_3['phis'][i]
        H_pred_4 = torsion_pred_data_4['phis'][i]
        H_pred_5 = torsion_pred_data_5['phis'][i]
        H_pred_6 = torsion_pred_data_6['phis'][i]
        H_pred_7 = torsion_pred_data_7['phis'][i]
        H_pred_8 = torsion_pred_data_8['phis'][i]
        
        r1 = torsion_relative_error(S_pred_1, H_pred_1, S, H)
        r2 = torsion_relative_error(S_pred_2, H_pred_2, S, H)
        r3 = torsion_relative_error(S_pred_3, H_pred_3, S, H)
        r4 = torsion_relative_error(S_pred_4, H_pred_4, S, H)
        r5 = torsion_relative_error(S_pred_5, H_pred_5, S, H)
        r6 = torsion_relative_error(S_pred_6, H_pred_6, S, H)
        r7 = torsion_relative_error(S_pred_7, H_pred_7, S, H)
        r8 = torsion_relative_error(S_pred_8, H_pred_8, S, H)
        #print(r1,r2,r3,r4,r5,r6,r7,r8)

        r = np.pi
        if r1 < r:
             r = r1
             torsion_pred_path = t1
             S_pred = S_pred_1
             H_pred = H_pred_1
        if r2 < r:
             r = r2
             torsion_pred_path = t2
             S_pred = S_pred_2
             H_pred = H_pred_2
        if r3 < r:
             r = r3
             torsion_pred_path = t3
             S_pred = S_pred_3
             H_pred = H_pred_3
        if r4 < r:
             r = r4
             torsion_pred_path = t4
             S_pred = S_pred_4
             H_pred = H_pred_4
        if r5 < r:
             r = r5
             torsion_pred_path = t5
             S_pred = S_pred_5
             H_pred = H_pred_5
        if r6 < r:
             r = r6
             torsion_pred_path = t6
             S_pred = S_pred_6
             H_pred = H_pred_6
        if r7 < r:
             r = r7
             torsion_pred_path = t7
             S_pred = S_pred_7
             H_pred = H_pred_7
        if r8 < r:
             r = r8
             torsion_pred_path = t8
             S_pred = S_pred_8 
             H_pred = H_pred_8
        print(p, torsion_pred_path, torsion_relative_error(S_pred, S, H_pred, H))

                       
import json

dir = "tests/DistanceMap/"

torsion_phi_start = [-1.5707963267948966,2.79252680319092] # always same values in gt
torsion_psi_end = [0.8726646259971648,4.1887902047863905]
torsion_omega_end = [0.017453292519943295,0.0]
def distance_and_torsion(pdbarg, distance_pred_path, torsion_pred_path, gt_path):
    
    # Change is_pred to True when passing in a valid prediction file!
    distance_pred_data = load_data(distance_pred_path, True, True, False)
    torsion_pred_data = load_data(torsion_pred_path, True, False, True)
    gt_data = load_data(gt_path, False, True, True)

    for (i, p) in enumerate(gt_data['pdbs']):
        # Torsion reconstruction.
        #tor_gt_coords = torsion_to_3d(gt_data['aas'][i], 
        #                gt_data['phis'][i], 
        #                gt_data['psis'][i])
        #tor_pred_coords = torsion_to_3d(gt_data['aas'][i], 
        #                pred_data['phis'][i], 
        #                pred_data['psis'][i])
        pred_dcalpha = distance_pred_data['dcalphas'][i] 
        pred_phis = torsion_pred_data['phis'][i]
        pred_psis = torsion_pred_data['psis'][i]
        
        gt_dcalpha = gt_data['dcalphas'][i]

        data = {}
        pdb = distance_pred_data['pdbs'][i]
        data['metaData'] = {'proteinName':pdb}
        data['distances'] = {}
        for j in range(0, len(pred_dcalpha)):
             prediction_j = pred_dcalpha[j]
             arr = [[float(prediction_j[k]), 0.0] for k in range(len(prediction_j))]
             #print(j, arr0)
             distance_j = {}
             for idx, val in enumerate(arr):
                 distance_j[str(idx)] = val
             #d0 = dict(enumerate(arr0))
             data['distances'][str(j)] = distance_j
#        print(data)
#        exit()
        data['proteinName'] = pdb
        aa = gt_data['aas'][i]
        data['proteinSequence'] = aa
        torsions = {}
        for j in range(0, len(pred_phis)):
             if pred_phis[j] == None:
                 pred_phis[j] = 0
             if pred_psis[j] == None:
                 pred_psis[j] = 0
             phi = pred_phis[j]*np.pi/180.0
             psi = pred_psis[j]*np.pi/180.0
             torsions[str(j)] = ([[phi,0.0],[psi,0.0],[np.pi,0.0]])
             if j==0:
                 torsions[str(j)] = ([torsion_phi_start,[psi,0.0],[np.pi,0.0]])
             if j==len(pred_phis)-1:
                 torsions[str(j)] = ([[phi,0.0],torsion_psi_end,torsion_omega_end])
        data['torsions'] = torsions
        if pdb == pdbarg:
            suffix = '.distance_torsions.json'
            filename = pdb+suffix
            with open(filename, 'w') as outfile:  
                json.dump(data, outfile, indent=1)
            reconstruct_torsion_and_distance(filename, pdb)

def reconstruct_torsion_and_distance(filename, pdb):
    print(pdb)
    pwd = '/Users/iddodrori/develop/protein/tsp/cuprotein/'
    path = '/Users/iddodrori/develop/protein/meshi9.36_light/out/tests/DistanceMap/'+pdb+'/'
    os.system("mv {} {}".format(filename, path))
    try:
        os.chdir(path)
    except:
        os.mkdir(path) 
        #os.system("cp {} {}".format('commands3', path)) 
        os.chdir(path)
    #proc = subprocess.Popen('java -cp ../../../ programs/DistanceMap2PDB {} commands3 refine 0'.format(filename), stdout=subprocess.PIPE, shell=True)
    #(out, err) = proc.communicate()
    os.chdir(pwd)
    return
 
def main2(pred_path, gt_path):
    # Change is_pred to True when passing in a valid prediction file!
    pred_data = load_data(pred_path, is_pred=True)
    gt_data = load_data(gt_path, is_pred=False)

    for (i, p) in enumerate(gt_data['pdbs']):
        # Torsion reconstruction.
        #tor_gt_coords = torsion_to_3d(gt_data['aas'][i], 
        #                gt_data['phis'][i], 
        #                gt_data['psis'][i])
        #tor_pred_coords = torsion_to_3d(gt_data['aas'][i], 
        #                pred_data['phis'][i], 
        #                pred_data['psis'][i])
        pred_dcalpha = pred_data['dcalphas'][i]
        gt_dcalpha = gt_data['dcalphas'][i]
        
        #tor_pred_data = gt_data['aas'][i], pred_data['phis'][i], pred_data['psis'][i]
        #gt_d = gt_data['aas'][i], gt_data['phis'][i], gt_data['psis'][i], gt_data['coords'][i]
        #print_metrics(p, gt_path, pred_path, gt_d, pred_dcalpha, method='dcalpha')
        """
        plt.imshow(gt_dcalpha, cmap='plasma')
        plt.savefig('testing_{}_dcalpha.png'.format(p))
        plt.imshow(pred_dcalpha, cmap='plasma')
        plt.savefig('{}_{}_dcalpha.png'.format(pred_path[:-4], p))
        """
        
        #tor_drmsd = calc_drmsd(tor_gt_coords, tor_pred_coords)
        #print_metrics(p, gt_path, pred_path, gt_d, tor_pred_data, method='torsion')

        #gt_coords = np.array(gt_data['coords'][i])
        #pred_coords = np.array(pred_data['coords'][i])
        #coord_drmsd = calc_drmsd(gt_coords, pred_coords)
        #coord_gdt_ts = calc_gdt_ts(gt_coords, pred_coords)

        #print("3dcoord, {}, {}, {}".format(p, coord_drmsd, coord_gdt_ts))
        #print("torsion, {}, {}, {}".format(p, tor_drmsd, tor_gdt_ts))

def get_amino_acid_sequences(pdb_path):
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path[:-4], pdb_path)
    ppb = PPBuilder()
    pdb_aas = []
    for pp in ppb.build_peptides(structure): 
        pdb_aa = str(pp.get_sequence())
        pdb_aas.append(pdb_aa)
    return pdb_aas
    
def pdb_to_torsion_and_res(pdb_path):
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path[:-4], pdb_path)
    A = []
    ppb = PPBuilder()
    pdb_aas = []
    model = ppb.build_peptides(structure)
    chain = model[0]
    phi_psi_list = chain.get_phi_psi_list()
    res_list = get_amino_acid_sequences(pdb_path)[0]
    return [x[0] for x in phi_psi_list], [x[1] for x in phi_psi_list], [x[0] for x in res_list]

def torsion_radians_to_degrees_sanitize(phi, psi):
    for i in range(0, len(phi)): # convert to degrees, first is None
        if phi[i] != None:
            phi[i] = phi[i] * (180/np.pi)
        else:
            phi[i] = 0
    for i in range(0, len(psi)): # convert to degrees, last is None
        if psi[i] != None:
            psi[i] = psi[i] * (180/np.pi)
        else:
            psi[i] = 0				
    return phi, psi

def torsion_and_res_to_pdb(phis, psis, ress, pdb_path):
	topo_filename = '{}.topo'.format(pdb_path)
	topo_file = open(topo_filename, 'w')
	for i in range(len(ress)):
		topo_file.write(str(phis[i])+' '+str(psis[i])+' '+str(ress[i])+'\n')
	topo_file.close()
	proc = subprocess.Popen('perl pdb_from_torsions.pl {} > {}'.format(topo_filename, pdb_path), stdout=subprocess.PIPE, shell=True)

def torsion_only_topo(pred_path, gt_path):
    pred_data = load_data(pred_path, is_pred=True)
    gt_data = load_data(gt_path, is_pred=False)
    for (i, p) in enumerate(gt_data['pdbs']):
        gt_filename = 'testing_{}.topo'.format(p)
        pred_filename = '{}_{}.topo'.format(pred_path[:-4], p)
        gt_file = open(gt_filename, 'w')
        pred_file = open(pred_filename, 'w')
        for j in range(len(gt_data['aas'][i])):
            gt_phi = gt_data['phis'][i][j]
            gt_psi = gt_data['psis'][i][j]
            res = gt_data['aas'][i][j]
            pred_phi = pred_data['phis'][i][j]
            pred_psi = pred_data['psis'][i][j]
            gt_file.write(str(gt_phi)+' '+str(gt_psi)+' '+str(res)+'\n')
            pred_file.write(str(pred_phi )+' '+str(pred_psi)+' '+str(res)+'\n')	
        gt_file.close()
        pred_file.close()
        gt_pdb_filename = 'testing_{}.pdb'.format(p)
        pred_pdb_filename = '{}_{}.pdb'.format(pred_path[:-4], p)
        proc = subprocess.Popen('perl pdb_from_torsions.pl {} > {}'.format(gt_filename, gt_pdb_filename), stdout=subprocess.PIPE, shell=True)
        proc = subprocess.Popen('perl pdb_from_torsions.pl {} > {}'.format(pred_filename, pred_pdb_filename), stdout=subprocess.PIPE, shell=True)

if __name__=='__main2__':
    parser = argparse.ArgumentParser()
    parser.add_argument("d1", help="Path to pickle file for distance predictions")
    parser.add_argument("d2", help="Path to pickle file for distance predictions")
    parser.add_argument("d3", help="Path to pickle file for distance predictions")
    parser.add_argument("d4", help="Path to pickle file for distance predictions")
    parser.add_argument("d5", help="Path to pickle file for distance predictions")
    parser.add_argument("d6", help="Path to pickle file for distance predictions")
    parser.add_argument("d7", help="Path to pickle file for distance predictions")
    parser.add_argument("d8", help="Path to pickle file for distance predictions")

    #parser.add_argument("t1", help="Path to pickle file for torsion predictions")
    parser.add_argument("gt", help="Path to pickle file for ground truth")
    args = parser.parse_args()
    d1 = args.d1
    d2 = args.d2
    d3 = args.d3
    d4 = args.d4
    d5 = args.d5
    d6 = args.d6
    d7 = args.d7
    d8 = args.d8

    #t1 = args.t1
    gt = args.gt
    
    # pdb -> torsion, res -> pdb
    #in_path = '2n64.pdb'
    #out_path = '2n64_recon.pdb'
    #phi, psi, res = pdb_to_torsion_and_res(in_path)
    #phi, psi = torsion_radians_to_degrees_sanitize(phi, psi)
    #torsion_and_res_to_pdb(phi, psi, res, out_path)
    #distance_and_torsion(d1, t1, gt)
    distance_only_comparison(d1, d2, d3, d4, d5, d6, d7, d8, gt)
    #distance_only(d1, gt)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    #parser.add_argument("pdb", help="pdb name")
    parser.add_argument("d1", help="Path to pickle file for distance predictions")
    #parser.add_argument("t1", help="Path to pickle file for torsion predictions")
    #parser.add_argument("t2", help="Path to pickle file for torsion predictions")
    #parser.add_argument("t3", help="Path to pickle file for torsion predictions")
    #parser.add_argument("t4", help="Path to pickle file for torsion predictions")
    #parser.add_argument("t5", help="Path to pickle file for torsion predictions")
    #parser.add_argument("t6", help="Path to pickle file for torsion predictions")
    #parser.add_argument("t7", help="Path to pickle file for torsion predictions")
    #parser.add_argument("t8", help="Path to pickle file for torsion predictions")

    parser.add_argument("gt", help="Path to pickle file for ground truth")
    args = parser.parse_args()
    #pdb = args.pdb
    d1 = args.d1
    #t1 = args.t1
    #t2 = args.t2
    #t3 = args.t3
    #t4 = args.t4
    #t5 = args.t5
    #t6 = args.t6
    #t7 = args.t7
    #t8 = args.t8

    gt = args.gt
    
    # pdb -> torsion, res -> pdb
    #in_path = '2n64.pdb'
    #out_path = '2n64_recon.pdb'
    #phi, psi, res = pdb_to_torsion_and_res(in_path)
    #phi, psi = torsion_radians_to_degrees_sanitize(phi, psi)
    #torsion_and_res_to_pdb(phi, psi, res, out_path)
    #distance_and_torsion(pdb, d1, t1, gt)
    #torsion_only_comparison(t1, t2, t3, t4, t5, t6, t7, t8, gt)
    distance_only(d1, gt)

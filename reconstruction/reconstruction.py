"""
	reconstruction.py
	
	Iddo Drori, Darshan Thaker
	Columbia University
	
	3D reconstructions from distance matrix
	
"""

import numpy as np
from scipy.spatial import distance_matrix
from sklearn import manifold
import cvxpy as cp

def mds(D):
	n = D.shape[0]
	dim = 3
	seed = np.random.RandomState(seed=3)
	mds = manifold.MDS(n_components = dim, max_iter = 3000, eps=1e-9, random_state=seed, dissimilarity="precomputed")
	Xtag = mds.fit(D).embedding_
	return Xtag

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
    prob = cp.Problem(cp.Minimize(obj), constraints)
    prob.solve(verbose=False)
    U, S, V = np.linalg.svd(G.value)
    u, s = np.real(U), np.real(S)
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
	prob.solve(verbose=False)
	U, S, V = np.linalg.svd(B.value)
	dim = 3
	Xtag = (np.sqrt(S)*np.transpose(V))[:,0:dim]
	return Xtag

def obj(x, args):
    n, D, Z, U, lam, rho = args
    G = x[:n*n].reshape(n, n)
    S = x[n*n:].reshape(n, n)
    F1 = lam * np.linalg.norm(S, 1)
    Gii = np.diag(G)
    F2 = Gii.reshape(n, 1) + Gii.reshape(1, n) - 2 * G + S - np.square(D)
    F2 = 0.5 * np.square(F2).sum()
    F3 = 0.5 * rho * (np.linalg.norm(G-Z+U, 2) ** 2)
    return F1 + F2 + F3
    
def gradient(x, args):
    n, D, Z, U, lam, rho = args
    G = x[:n*n].reshape(n, n)
    S = x[n*n:].reshape(n, n)
    Gii = np.diag(G)
    dLOSS = Gii.reshape(n,1) + Gii.reshape(1,n) - 2*G + S - D**2
    dG = np.copy(dLOSS)
    np.fill_diagonal(dG, 0)
    dGii = dG.sum(0) + dG.sum(1)    
    dG = -2*dG + np.diag(dGii)
    np.fill_diagonal(dG, 0)
    dG += rho*(G-Z+U)    
    dS = lam * np.sign(S) + dLOSS
    return np.r_[dG, dS].ravel()

def gradient_descent_with_momentum(x, args, alpha, beta):
    it, step, max_iters = 0, 10, 100
    m = np.zeros_like(x)
    while step > 0.1 and it < max_iters:
        it += 1
        x_ = np.copy(x)
        g = gradient(x, args)
        m = beta * m + alpha * g
        x = x - m
        step = np.abs(obj(x, args) - obj(x_, args))
    return x

def gram_matrix(x, args):
    it, step = 0, 100
    alpha, beta = 0.03, 0.95
    while step > 0.1:
        it += 1
        x_ = np.copy(x)
        n, D, Z, U, lam, rho = args
        x = gradient_descent_with_momentum(x, args, alpha, beta)
        G = x[:n*n].reshape(n, n)
        S = x[n*n:].reshape(n, n)
        l, v = np.linalg.eig(G + U)
        idx = l.argsort()[::-1]
        l = l[idx]
        l = (l > 0) * l
        v = v[:, idx]
        Z = np.dot(np.dot(v, np.diag(l)), v.T)
        U = U + G - Z
        F_ = obj(x_, args)
        args = (n, D, Z, U, lam, rho)
        F = obj(x, args)
        step = np.abs(F - F_)
    return G

def admm(D):
    lam = 1
    rho = 100
    n = D.shape[0]
    Z = np.ones_like(D)
    U = np.ones_like(D)
    G = np.ones_like(D)
    S = np.ones_like(D)
    x = np.r_[G, S].ravel()
    args = (n, D, Z, U, lam, rho)
    G = gram_matrix(x, args)
    U, S, V = np.linalg.svd(G)
    u, s = np.real(U), np.real(S)
    dim = 3
    Xtag = np.dot(u, np.diag(np.sqrt(s)))[:,0:dim]
    return Xtag

def relative_error(A, Atag):
    return np.linalg.norm(Atag - A)/np.linalg.norm(A)

n = 30
dim = 3
X = np.random.randn(n, dim)
D = distance_matrix(X, X)
Xtag = mds(D)
print('mds relative error', relative_error(D, distance_matrix(Xtag, Xtag)))
Xtag = sdp(D)
print('sdp relative error', relative_error(D, distance_matrix(Xtag, Xtag)))
Xtag = sdr(D)
print('sdr relative error', relative_error(D, distance_matrix(Xtag, Xtag)))
Xtag = admm(D)
print('admm relative error', relative_error(D, distance_matrix(Xtag, Xtag)))



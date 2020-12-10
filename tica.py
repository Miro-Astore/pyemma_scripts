#import pyemma.coordinates
import pyemma
import numpy as np
from pylab import *
import pdb
matplotlib.rcParams.update({'font.size': 14})


import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
# Switch progress bars off, because this notebook would generate a lot of progress bars
pyemma.config.show_progress_bars = 'True'


def project_and_cluster(trajfiles, featurizer, sparsify=False, tica=True, lag=100000, scale=True, var_cutoff=1.0 , ncluster=100):
    """
    Returns
    -------
    trans_obj, Y, clustering

    """
    X = coor.load(trajfiles, featurizer)
    if sparsify:
        X = remove_constant(X)
    if tica:
        trans_obj = coor.tica(X, lag=lag, var_cutoff=var_cutoff)
        Y = trans_obj.get_output()
    else:
        trans_obj = coor.pca(X, dim=-1, var_cutoff=var_cutoff)
        Y = trans_obj.get_output()
    if scale:
        for y in Y:
            y *= trans_obj.eigenvalues[:trans_obj.dimension()]
    if cluster:
        cl_obj = coor.cluster_kmeans(Y, k=ncluster, max_iter=3, fixed_seed=True)
        return trans_obj, Y, cl_obj
    return trans_obj, Y

def eval_transformer(trans_obj):
    # Effective dimension (Really? If we just underestimate the Eigenvalues this value also shrinks...))
    print('Evaluating transformer: ', str(trans_obj.__class__))
    print('effective dimension', np.sum(1.0 - trans_obj.cumvar))
    print('eigenvalues', trans_obj.eigenvalues[:5])
    print('partial eigensum', np.sum(trans_obj.eigenvalues[:10]))
    print('total variance', np.sum(trans_obj.eigenvalues ** 2))
    print()

def plot_map(Y, sx=None, sy=None, tickspacing1=1.0, tickspacing2=1.0, timestep=1.0, timeunit='ns'):
    if not isinstance(Y, np.ndarray):
        Y = Y[0]
    if sx is None:
        sx = -np.sign(Y[0,0])
    if sy is None:
        sy = -np.sign(Y[0,1])
    Y1 = sx*Y[:, 0]
    min1 = np.min(Y1)
    max1 = np.max(Y1)
    Y2 = sy*Y[:, 1]
    min2 = np.min(Y2)
    max2 = np.max(Y2)
    # figure
    figure(figsize=(16,4))
    # trajectories
    subplot2grid((2,2), (0,0))
    plot(timestep*np.arange(len(Y1)), Y1)
    xlim(0, timestep*len(Y1))
    yticks(np.arange(int(min1), int(max1)+1, tickspacing1))
    ylabel('component 1')
    subplot2grid((2,2), (1,0))
    plot(timestep*np.arange(len(Y2)), Y2)
    xlim(0, timestep*len(Y2))
    ylabel('component 2')
    yticks(np.arange(int(min2), int(max2)+1, tickspacing2))
    xlabel('time / ' + timeunit)
    # histogram data
    subplot2grid((2,2), (0,1), rowspan=2)
    z,x,y = np.histogram2d(Y1, Y2, bins=50)
    z += 0.1
    # compute free energies
    F = -np.log(z)
    # contour plot
    extent = [x[0], x[-1], y[0], y[-1]]
    xticks(np.arange(int(min1), int(max1)+1, tickspacing1))
    yticks(np.arange(int(min2), int(max2)+1, tickspacing2))
    contourf(F.T, 50, cmap=plt.cm.nipy_spectral, extent=extent)
    xlabel('component 1')
    ylabel('component 2')

#feat=pyemma.coordinates.featurizer('prot.pdb')


top = '/scratch/f91/ma2374/vsite_CFTR/wt/310K/combined_pca_analysis/cov-domain-average.pdb'
#trajs = ['/scratch/f91/ma2374/vsite_CFTR/wt/310K/1/wt_ca.xtc','/scratch/f91/ma2374/vsite_CFTR/wt/310K/2/wt_ca.xtc','/scratch/f91/ma2374/vsite_CFTR/wt/310K/3/wt_ca.xtc','/scratch/f91/ma2374/vsite_CFTR/wt/310K/combined_pca_analysis/aa_wt_ca.xtc']
trajs = ['/scratch/f91/ma2374/vsite_CFTR/wt/310K/1/wt_ca_domain.xtc','/scratch/f91/ma2374/vsite_CFTR/wt/310K/2/wt_ca_domain.xtc','/scratch/f91/ma2374/vsite_CFTR/wt/310K/3/wt_ca_domain.xtc','/scratch/f91/ma2374/vsite_CFTR/wt/310K/combined_pca_analysis/aa_wt_ca_domain.xtc']

feat_Ca = coor.featurizer(top)
feat_Ca.add_selection(feat_Ca.select ('name CA'))
print(feat_Ca.dimension())
cluster = False
if cluster:
    tica_Ca, tica_Y_Ca, tica_cl_Ca = project_and_cluster(trajs, feat_Ca)
else:
    tica_Ca, tica_Y_Ca = project_and_cluster(trajs, feat_Ca)
    #tica_Ca, tica_Y_Ca = project_and_cluster(trajs, feat_Ca)
print(np.shape(tica_Ca.eigenvectors))
x=(tica_Ca.get_params())
#pdb.set_trace()
np.save('tica_eigvec.npy',tica_Ca.eigenvectors)
np.save('tica_eigval.npy',tica_Ca.eigenvalues)
print('feat_means.npy',tica_Ca.get_params().keys())
#pdb.set_trace()
#np.save('feat_means.npy',tica_Ca.get_params()['mean'])
#eval_transformer(tica_Ca)
#pdb.set_trace()
#pca_Ca, pca_Y_Ca, pca_cl_Ca = project_and_cluster(trajs, feat_Ca, tica=False)
#eval_transformer(pca_Ca)

feat_CA.add_distances_ca()
reader=pyemma.coordinates.source(['prot.xtc'],features=feat)
tic = pyemma.coordinates.tica(reader,lag=2,dim=1,chunksize=1,ncov_max=1).get_output()[0]
print (tic)
np.savetxt('tica.dat',tic)
#

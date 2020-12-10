import numpy as np 
import MDAnalysis as mda

#import pdb and construct extract cas.
u=mda.Universe('prot.pdb')
cas=u.select_atoms('name CA')
#tica_corr=np.load('./tica_eigvec.npy')
tics=np.load('./tica_eigvec.npy')
eigvals=np.load('./tica_eigval.npy')
feat_means=np.load('./feat_means.npy')
#get positions and reshape into 1*3n array so we can use the tica corr matrix to act on it.
positions_list=np.array(cas.positions)
#tot_elements=np.shape(positions_list)[0]*np.shape(positions_list)[1]
#pos_vec=np.reshape(positions_list,(tot_elements,1))
#perform transformation
#new_coors_vec = (tica_corr.dot(pos_vec))
#why does this work? nobody knows. eigen vectors are dmeaned so we have to reintroduce the mean. units are in nm so convert them to AA for vmd
for i in range(5):
    new_coors_vec = (tics[:,i]+feat_means.T)*10
    new_coors_list = (np.reshape(new_coors_vec,(np.shape(positions_list)[0],3)))
    cas.positions=new_coors_list
    #cas.write('moment_of_truth.pdb')
    cas.write((str (i) + '.pdb'))

    mean_coors_vec = (feat_means.T)*10
    mean_coors_list = (np.reshape(mean_coors_vec,(np.shape(positions_list)[0],3)))
cas.positions=mean_coors_list
cas.write('start.pdb')

import numpy as np 
import MDAnalysis as mda
<<<<<<< HEAD

#import pdb and construct extract cas.
u=mda.Universe('prot.pdb')
cas=u.select_atoms('name CA')
=======
import pdb

#import pdb and construct extract cas.
u=mda.Universe('bpti_prot.pdb')
cas=u.select_atoms('name CA')
#pdb.set_trace()
>>>>>>> cb5bbc6ba4c2aa07ae080cfca01cd1d52466770e
#tica_corr=np.load('./tica_eigvec.npy')
tics=np.load('./tica_eigvec.npy')
eigvals=np.load('./tica_eigval.npy')
feat_means=np.load('./feat_means.npy')
<<<<<<< HEAD
#get positions and reshape into 1*3n array so we can use the tica corr matrix to act on it.
positions_list=np.array(cas.positions)
=======
#feat_means=np.load('./feat_means.npy')
#get positions and reshape into 1*3n array so we can use the tica corr matrix to act on it.
#positions_list=(np.matrix(cas.positions))
#positions_list=np.array(positions_list.flatten())  
#positions_list=positions_list[0]

positions_list = feat_means[0]
>>>>>>> cb5bbc6ba4c2aa07ae080cfca01cd1d52466770e
#tot_elements=np.shape(positions_list)[0]*np.shape(positions_list)[1]
#pos_vec=np.reshape(positions_list,(tot_elements,1))
#perform transformation
#new_coors_vec = (tica_corr.dot(pos_vec))
#why does this work? nobody knows. eigen vectors are dmeaned so we have to reintroduce the mean. units are in nm so convert them to AA for vmd
for i in range(5):
<<<<<<< HEAD
    new_coors_vec = (tics[:,i]+feat_means.T)*10
    new_coors_list = (np.reshape(new_coors_vec,(np.shape(positions_list)[0],3)))
=======
    #new_coors_vec = (tics[:,i]/(np.linalg.norm(tics[:,i])*eigvals[i]))*10+feat_means.T
    inv_tic_vec=np.linalg.inv(tics)
    basis_vec=np.zeros(len(eigvals))
    basis_vec[0] = 1
    add_vec = basis_vec.dot(inv_tic_vec)
    new_coors_vec = 10*((add_vec)*2+feat_means.T)
    #new_coors_vec = (tics[:,i])*10
    print(positions_list[0])
    new_coors_list = (np.reshape(new_coors_vec,(((int(len(positions_list)/3),3)))))
>>>>>>> cb5bbc6ba4c2aa07ae080cfca01cd1d52466770e
    cas.positions=new_coors_list
    #cas.write('moment_of_truth.pdb')
    cas.write((str (i) + '.pdb'))

<<<<<<< HEAD
    mean_coors_vec = (feat_means.T)*10
    mean_coors_list = (np.reshape(mean_coors_vec,(np.shape(positions_list)[0],3)))
=======
    mean_coors_vec = (feat_means.T)
    mean_coors_list = (np.reshape(len(cas),(np.shape(positions_list)[0],3)))
    #mean_coors_list = (np.reshape(mean_coors_vec,(((int(len(positions_list)/3),3)))))
>>>>>>> cb5bbc6ba4c2aa07ae080cfca01cd1d52466770e
cas.positions=mean_coors_list
cas.write('start.pdb')

import pdb
import numpy as np
import MDAnalysis as mda
u = mda.Universe('prot.psf', 'out.xtc')

# selection of atomgroups
ags = [res.phi_selection() for res in u.residues[4:9]]

from MDAnalysis.analysis.dihedrals import Dihedral
from MDAnalysis.analysis.dihedrals import Ramachandran
#R = Dihedral(ags).run()
r=u.select_atoms("resid 5-10")
R = Ramachandran(r).run()

#R = Dihedral(ags).run()

np.save('out2.npy',R.angles)

#R.save()

#import matplotlib.pyplot as plt

#fig, ax = plt.subplots(figsize=plt.figaspect(1))

#R.plot(ax=ax, color='k', marker='s')
#pdb.set_trace()

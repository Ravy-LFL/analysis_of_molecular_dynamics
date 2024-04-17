"""
    This script aim to compute distance between residue along a trajectory.
    The results are write in a csv file.
"""

import numpy as np
from MDAnalysis.analysis import distances
import MDAnalysis as mda
import sys

__author__ = "Ravy LEON FOUN LIN"
__version__ = "1.0.0"
__date__ = "15/02/2024"

STRUCTURE = sys.argv[1]

#  Path of the needed files.
PATH = STUCTURE

#  Coordinate file.
GRO = PATH+"protein.pdb"

#  Trajectory file.
XTC = PATH+"fit.xtc"


def compute_distance(u, COUNT_RESIDUE, frame,file) :
    """
        Compute distance between each residue (one vs all).
        Write the csv file which contain the distances.

        Parameters
        ----------
            u = mda.Universe()
            COUNT_RESIDUE = Number of residues for one protein.
            frame = number of the frame.
            file = csv file.

        Returns
        -------
            int : 0
    """

    #  Iterate on the residue number for the protein A.
    for i in range(1,COUNT_RESIDUE) :
        sel_A = f'name BB and resid {i} and segid A'

        #  Iterate on the residue number, starting by i, for the protein B.
        for j in range(i,COUNT_RESIDUE) :

            #  We skip the residue if it is neighbour to the actual residue in A.
            if j == i+1 :
                continue
            else :
                #  Select the backbone grain of protein B.
                sel_B = f'name BB and resid {j} and segid B'

                #  Select the residues.
                residue_i = u.select_atoms(sel_A)
                residue_j = u.select_atoms(sel_B)
                
                #  Compute the distances.
                dist = distances.dist(residue_i,residue_j)

                #  Retrieve the variables.
                id_i = str(dist[0][0])
                id_j = str(dist[1][0])
                name_i = residue_i.resnames[0]
                name_j = residue_j.resnames[0]
                dist_val = str(dist[2][0])

                #  Write the data.
                file.write(f"{id_i},{id_j},{name_i},{name_j},{dist_val},{frame}\n")

    return 0


if __name__ == '__main__' :

    #  Define universe object.
    u = mda.Universe(GRO,XTC)

    #  Count of residue in one protein.
    COUNT_RESIDUE = int(u.residues.n_residues/2)

    #  Open the csv file.
    file = open(f'./{STRUCTURE}_count.csv','w')
    file.write("resi_i,resi_j,name_i,name_j,distance,frame\n")

    #  Frame number.
    frame = 1

    #  Iterate on the trajectory.
    print('count...')
    for ts in u.trajectory :

        #  Compute data.
        compute_distance(u, COUNT_RESIDUE, frame,file)
    
        #  Increment frame count.
        frame +=1

    #  Close the file.
    file.close()

    print('all done')

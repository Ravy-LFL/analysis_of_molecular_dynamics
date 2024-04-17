""" Generate X, Y and labels value from the feather files generated grom trajectories.
    
    Usage
    =====
        python3 script_generate.py -file <path>

"""

#  Import libraries.
import argparse
import pandas as pd
import pickle as pkl
import sys

#  Informations.
__author__ = "Ravy LEON FOUN LIN"

#  Parsing.
parser = argparse.ArgumentParser()
parser.add_argument("-file", help="Path to the feather file")
args = parser.parse_args()
FILE = args.file

#  Check if parsing is okay.
if FILE == None :
    parser.print_help()
    sys.exit()



def build_dict_resi_frame(interaction_df, r_interactif):
    """Build a dictionnary which have as key the residues and as value a list of frame.

        Parameters
        ----------
            interaction_df (pandas.DataFrame) : Contains the residues which interact in the simulation with other informations.
            r_interactif (list) : List of residues (ID) which interact.

        Returns
        -------
            Dictionnary of residue in interactions and the frame where they interact.
    """

    #  Dictionnary return.
    dict_resi_frame = dict()


    #  Iterate on the DataFrame of residues interacting.
    for index, row in interaction_df.iterrows() :
        
        #  Retrieve ID of the first residue.
        resi_i = row['resi_i']

        #  Retrieve ID of the second residue.
        resi_j = row['resi_j']

        #  Retrieve the frame.
        frame = row['frame']
        
        #  Check if residue is in the list of residue in interaction.
        if resi_i in r_interactif :
            
            # Add it to the dictionnary as a key + the frame associated.
            if resi_i not in dict_resi_frame.keys() :
                dict_resi_frame[resi_i] = list()
                dict_resi_frame[resi_i].append(frame)
            
            #  If already in dict, just add the frame to the existing list.
            else :
                dict_resi_frame[resi_i].append(frame)
                
        #  Same treatment than for the first residue.
        if resi_j in r_interactif :
        
            if resi_j not in dict_resi_frame.keys() :
                dict_resi_frame[resi_j] = list()
                dict_resi_frame[resi_j].append(frame)
            
            else :
                dict_resi_frame[resi_j].append(frame)
    
    #  Return the dictionnary.
    return dict_resi_frame

def generate_xy(path : str) :
    """Generate the x, y and labels value for the feather file representing a trajectory.

        Parameters
        ----------
            path (str) : Path to the feather file.
        
        Returns
        -------
            float & str : List of index and the associated labels.
    """

    #  Y values (the frames).
    y = list()

    #  X values (the residue IDs).
    x = list()

    #  Labels name (residue names).
    labels = list()

    #  Create the dataframe from the feather.
    df_iter = pd.read_csv(path, sep = ',', chunksize=50000)
    
    for df in df_iter : 
        #  Only retrieve the ID and the name associated for each residues.
        df_resi_type = df[['resi_i','name_i']].drop_duplicates()

        #  Retrieve the pair at a distance lower than 5 A.
        interaction_df = df.query("distance <= 7")

        #  Create list of interacting residues, retrieve first residue ID of the pair.
        r_interacting = list(interaction_df['resi_i'].unique())

        #  Add to the list the second residue ID of the pair.
        r_interacting += [i for i in interaction_df['resi_j'].unique()]

        #  Create a dict which associate the residues to the frames.
        dict_resi_frame = build_dict_resi_frame(interaction_df, r_interacting)

        #  Iterate on the dictionnary of residue-frames association.
        for i in list(dict_resi_frame.keys()) :
            
            try :
                #  Retrieve the name of the residue in interaction.
                l = list(df_resi_type.query(f"resi_i == {i}")['name_i'])[0]
            except :
                continue
            #  Add them to the lists.
            for j in dict_resi_frame[i] :

                labels.append(l)
                x.append(i)
                y.append(j)

    return (x,y,labels)

if __name__ == '__main__' :

    #  Generate the values.
    print('run...')
    x,y,labels = generate_xy(FILE)

    NAME = FILE.split('/')[-1].split('_')[0]
    
    data = {'x':x, 'y':y, 'labels':labels}
    
    with open(f"{NAME}_data.pkl", 'wb') as f :
        pkl.dump(data,f)
    

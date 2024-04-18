# scripts_contacts_cg_martini3

### `compute_contact_along_trj.py`
This python script is use to generate the csv file which contains informations on the distances between the backbone beads along the trajectory.
Usage : `python3 compute_contact_along_trj.py <path_to_the_dir_data>`
The csv file will be write in this directory.

### `generate_axes_chunked.py`
This script use the csv file previously made and will create a pickle file with the residue ID, the frames where the residue is in a contact and the residue name.
It keeps every residue in a contact of at least 7 Angstrom.

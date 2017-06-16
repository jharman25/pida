__description__ = \
"""
A class for mapping protein-protein interaction hotspots onto crystal
structures based on results from docking programs (Autodock Vina,
ClusPro, Haddock, Swarmdock, Pie-Dock, Rosetta/PyRosetta, etc). 
"""

__author__ = "Joseph Harman"
__date__ = "2017-06-02"

from matplotlib import pyplot as plt
import numpy as np
import mdtraj as md
import nglview as nv
import os, itertools, pickle

class pida(object):
    """Class for parsing and visualizing Protein-protein Interaction Docking Analysis data. 
    
    
    Attributes:
        path: directory to data files as a string.
        new_path: new folder (string) gets made in path that will house files converted to mdtraj format.
        receptor_chains: list of integers corresponding to receptor chains. 
        ligand_chains: list of integers corresponding to ligand chains.
    """
    
    def __init__(self, path, new_path, receptor_chains, ligand_chains):
        
        self.path = path
        self.new_path = new_path
        self.ligand_chains = ligand_chains
        self.receptor_chains = receptor_chains
       
        self.mdtraj_converter = self.mdtraj_converter()
        self.aggregate_contacts, self.aggregate_residues = self.contacts_iterator()
        
        self.receptor_residues, self.receptor_counts, self.ligand_residues, self.ligand_counts, self.actual_ligand_residues = self.docking_histograms()
        
        self.receptor_atoms, self.receptor_atom_counts, self.actual_ligand_atoms, self.ligand_atoms, self.ligand_atom_counts = self.atom_contacts_parser()
        
        self.heatmap_viewer = self.heatmap()
        
        
    def mdtraj_converter(self, filetype = "ClusPro"):
        """
        Iterates through files in a user-specified directory and formats them for mdtraj input, 
        allowing mdtraj to read the files in as single-frame trajectories. Writes out new set of files
        labeled modeln.pdb (n = 0...n)
        """
        files = os.listdir(self.path)
        
        if (self.new_path) not in files:
            os.makedirs(self.path + self.new_path)
        
        new_titles = []

        for i in range(len(files)):
            i = "model%s.pdb" % (i) 
            new_titles.append(i)

        for i in range(len(files)):
            
            if filetype == "ClusPro":
                
                def cluspro_converter(self, filename, output):
                    """
                    Re-formats cluspro pdb files so mdtraj can read them in as single-frame trajectories.  
                    """
                
                    f = open(filename)
                    g = open(output, "w")

                    for line in f.readlines():
                        if line[:3] != "END":
                            g.write(line)

                    f.close()
                    g.close()
                    
                cluspro_converter(self, self.path + files[i], self.path + self.new_path + new_titles[i])
                
        print("Files successfully converted...")
        

    def contacts_iterator(self):
        """
        Iterates through all files in a user-specified path (that have been edited by function mdtraj_converter 
        for mdtraj compatibility), calculates interface contacts using mdtraj_compute_contacts 
        (here called by function get_interface_contacts), pulls out all atom indices in each interface 
        and adds all interface atoms and interface residues to aggregate list of contacts. 
        Aggregate_contacts gets printed out into a file named "aggregate_contacts.txt" for later use, if needed. 
        Not currently mapped directly to file names, but aggregate_contacts[i] corresponds to all atom contacts within model[i].pdb.
        aggregate_contacts[i] feeds directly into model_viewer for viewing of individual docking solutions, with docking
        interface highlighted.    
        """
        
        files = os.listdir(self.path + self.new_path)
        
        aggregate_contacts = []
        aggregate_residues = []
        
        def get_interface_contacts(frame, ca_cutoff_ang=10.):
            """
            Identify interface residues between ligand chains and receptor chains using mdtraj. 
            Residues identified by user-specified c-alpha cutoff, preset to 10 angstroms. 
            Feeds into contacts_iterator. 
            """

            #Get list of residues in receptor and ligand
            r_residues = []
            for chain in self.receptor_chains:
                r_residues.extend([residue.index for residue in frame.topology.chain(chain).residues])

            l_residues = []
            for chain in self.ligand_chains:
                l_residues.extend([residue.index for residue in frame.topology.chain(chain).residues])

            # Make an array of potential contact pairs between receptor and ligand
            contact_pairs = np.array([(i,j) for i in r_residues for j in l_residues])

            # Check which ones fall within c-alpha distance cutoff    
            is_contact = (10.*md.compute_contacts(frame, scheme='ca', contacts=contact_pairs)[0] < ca_cutoff_ang)[0]

            # Go from bool truth values to the actual residues
            contacts = contact_pairs[is_contact]

            # Go from pairs to flattened list of unique residues involved in contacts
            self.interface_residues = sorted(set(contacts[:,0]).union(set(contacts[:,1])))

            return self.interface_residues

        for i in range(len(files)):
            i = "model%s.pdb" % (i) 
            pdb = md.load(self.path  + self.new_path + i)
            res = get_interface_contacts(pdb)
            aggregate_residues.append(res)

            atom_indices = []

            for i in res:
                a = [atom.index for atom in pdb.topology.residue(i).atoms]
                atom_indices.append(a)

            interface_atoms = np.hstack(atom_indices).tolist()
            aggregate_contacts.append(interface_atoms)

        with open(self.path + self.new_path + "aggregate_contacts.pickle", "wb") as g:
            pickle.dump(aggregate_contacts, g)

        g.close()

        with open(self.path + self.new_path + "aggregate_residues.pickle", "wb") as h:
            pickle.dump(aggregate_residues, h)

        h.close()
        
        print("Contacts iterated...")

        return aggregate_contacts, aggregate_residues  

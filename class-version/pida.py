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
        self.view_single_solution = self.model_viewer()
        
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


    def model_viewer(self, model_number = 0):
        """
        Quick viewer for visually inspecting docking contacts output. 
        Takes aggregate_contacts output from contacts_iterator and plots contacts (in red) 
        onto receptor-ligand docking solution (blue). 
        """

        pdb = md.load(self.path  + self.new_path + "model%s.pdb" % (model_number))

        view_single_solution = nv.show_mdtraj(pdb)
        view_single_solution.add_cartoon(selection ="protein", color = "blue")
        view_single_solution.add_ball_and_stick(selection = self.aggregate_contacts[model_number], color = "red")
        
        print("Single models ready for viewing...")

        return view_single_solution


    def docking_histograms(self):
        """
        Takes output aggregate_residues from contacts_iterator and plots histograms of docking hits by residue
        for both receptor and ligand. Requires path name for directory in which all docking pdb files are stored,
        as well as the list of l_chains you used for contacts_iterator - this allows the function to    differentiate 
        between ligand and receptor residues. 
        """

        counts = {}

        for i in itertools.chain(*self.aggregate_residues):
            counts[i] = counts.get(i, 0) + 1

        all_resis = list(counts.keys())
        all_counts = list(counts.values())

        pdb = md.load(self.path + self.new_path + "model0.pdb")

        ligand_residue_1 = [residue.index for residue in pdb.topology.chain(self.ligand_chains[0]).residues][0]

        receptor_residues = []
        ligand_residues = []
        receptor_counts = []
        ligand_counts = []
        actual_ligand_residues = []

        for i in range(len(all_resis)):
            if all_resis[i] < ligand_residue_1:
                receptor_residues.append(all_resis[i])
                receptor_counts.append(all_counts[i])
            else:
                ligand_residues.append(all_resis[i])
                ligand_counts.append(all_counts[i])

        for i in range(len(ligand_residues)):
            actual_ligand_residues.append(ligand_residues[i]-ligand_residue_1)

        rec_fig = plt.figure()
        ax1 = rec_fig.add_subplot(1,1,1)
        ax1.set_xlabel("Receptor residue #")
        ax1.set_ylabel("Counts")
        ax1.bar(receptor_residues, receptor_counts)
        

        lig_fig = plt.figure()
        ax2 = lig_fig.add_subplot(1,1,1)
        ax2.set_xlabel("Ligand residue #")
        ax2.set_ylabel("Counts")
        ax2.bar(actual_ligand_residues, ligand_counts, color ="g")      
        
        print("Histograms constructed...")

        return receptor_residues, receptor_counts, ligand_residues, ligand_counts, actual_ligand_residues


    def atom_contacts_parser(self):
        """
        Parses through aggregate_contacts data from contacts_iterator and spits out 
        contacted receptor_atoms and ligand_atoms, and their corresponding counts. 
        This data feeds into the heatmap viewer. 
        """

        counts = {}

        for i in itertools.chain(*self.aggregate_contacts):
                counts[i] = counts.get(i, 0) + 1

        all_contacts = list(counts.keys())
        all_counts = list(counts.values())

        pdb = md.load(self.path + self.new_path + "model0.pdb")

        ligand_atom_1 = [atom.index for atom in pdb.topology.chain(self.ligand_chains[0]).atoms][0]

        receptor_atoms = []
        ligand_atoms = []
        actual_ligand_atoms = []
        receptor_atom_counts = []
        ligand_atom_counts = []

        for i in range(len(all_contacts)):
            if all_contacts[i] < ligand_atom_1:
                receptor_atoms.append(all_contacts[i])
                receptor_atom_counts.append(all_counts[i])
            else:
                ligand_atoms.append(all_contacts[i])
                ligand_atom_counts.append(all_counts[i])

        receptor_atoms = np.hstack(receptor_atoms).tolist()
        ligand_atoms = np.hstack(ligand_atoms).tolist()
        receptor_atom_counts = np.hstack(receptor_atom_counts).tolist()
        ligand_atom_counts = np.hstack(ligand_atom_counts).tolist()

        for i in range(len(ligand_atoms)):
            actual_ligand_atoms.append(ligand_atoms[i]-ligand_atom_1)

        actual_ligand_atoms = np.hstack(actual_ligand_atoms).tolist()  
        
        print("atom contacts parsed...")

        return receptor_atoms, receptor_atom_counts, actual_ligand_atoms, ligand_atoms, ligand_atom_counts

    
    def heatmap(self, plot_receptor = True, plot_ligand = False):
        """
        Constructs a heatmap of contact points onto receptor and ligand structures.
        Requires input from atom_contacts_parser. Color cutoffs are currently hard-coded
        for this particular data set, but can be easily modified. 
        """
        pdb = md.load(self.path + self.new_path + "model0.pdb")

        rec = []

        for i in range(len(self.receptor_chains)):
            c = [atom.index for atom in pdb.topology.chain(self.receptor_chains[i]).atoms]
            rec.append(c)

        rec = np.array(rec)
        rec = np.hstack(rec).tolist()
        rec_structure = md.load(self.path + self.new_path + "model0.pdb", atom_indices=rec)

        lig = []

        for i in range(len(self.ligand_chains)):
            c = [atom.index for atom in pdb.topology.chain(self.ligand_chains[i]).atoms]
            lig.append(c)

        lig = np.array(lig)
        lig = np.hstack(lig).tolist()
        lig_structure = md.load(self.path + self.new_path + "model0.pdb", atom_indices=lig)

        if plot_receptor == True:

            heatmap_viewer = nv.show_mdtraj(rec_structure)
            heatmap_viewer.clear()
            heatmap_viewer.add_ball_and_stick(selection ="residue", color = "blue")

            low = []
            med = []
            hi = []

            for i in range(len(self.receptor_atoms)):
                if self.receptor_atom_counts[i] > 1 and self.receptor_atom_counts[i] < 6:
                    low.append(self.receptor_atoms[i])
                elif self.receptor_atom_counts[i] > 5 and self.receptor_atom_counts[i] < 11:
                    med.append(self.receptor_atoms[i])
                else:
                    hi.append(self.receptor_atoms[i])

            heatmap_viewer.add_ball_and_stick(selection = low, color = "yellow")
            heatmap_viewer.add_ball_and_stick(selection = med, color = "orange")
            heatmap_viewer.add_ball_and_stick(selection = hi, color = "red")

        if plot_ligand == True:

            heatmap_viewer = nv.show_mdtraj(lig_structure)
            heatmap_viewer.clear()
            heatmap_viewer.add_ball_and_stick(selection ="residue", color = "blue")

            l = []
            m = []
            h = []

            for i in range(len(self.actual_ligand_atoms)):
                if self.ligand_atom_counts[i] > 6 and self.ligand_atom_counts[i] < 26:
                    l.append(self.actual_ligand_atoms[i])
                elif self.ligand_atom_counts[i] > 25 and self.ligand_atom_counts[i] < 51:
                    m.append(self.actual_ligand_atoms[i])
                else:
                    h.append(self.actual_ligand_atoms[i])

            heatmap_viewer.add_ball_and_stick(selection = l, color = "yellow")
            heatmap_viewer.add_ball_and_stick(selection = m, color = "orange")
            heatmap_viewer.add_ball_and_stick(selection = h, color = "red")
            
        print("Heatmaps constructed...done")

        return heatmap_viewer
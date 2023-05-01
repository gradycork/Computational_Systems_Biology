#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 17:35:44 2023

@author: gradycorkum
"""

# Import cell (RUN THIS CELL TO TEST YOUR INSTALL/ENVIRONMENT)
# part 1
from Bio.PDB import *
import nglview as nv
import ipywidgets
# part 2
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

#%% Show 3D model
pdb_parser = PDBParser()
structure = pdb_parser.get_structure("P", "P7_Structure.pdb")
view = nv.show_biopython(structure)
view
# TI_Thesis
Code and data for Tinbergen Institute's MPhil thesis
This is quick overview of the code files required for running the programs. They do not include the required dataset, but these are available upon request (or once I've got them uploaded in Erasmus archives).

The main file is TI_Thesis_Draft.py. This includes the bulk of the code and requires grad.py, AggregateACS_DRAFT.py, and AdjacencyMatrix_DRAFT.py. 

Expendnm.prg, expendnm.run, auxpanel.dat, and panel.dat are required for GAUSS estimation. The code for specific models that can be plugged into expendnm.prg can be found in gmm_x.text and lfm_x.txt. To swap panel.dat to other .csv files with the same structure, use the code in load_in_routine.txt.

Neighbour_working.py is code to extract neighbour relationships from GIS. The code was tested in QGIS 3.40.4.

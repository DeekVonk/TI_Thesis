# TI_Thesis
Code and data for Tinbergen Institute's MPhil thesis
This is quick overview of the code files required for running the programs. They do not include the required dataset, but these are available upon request (or once I've got them uploaded in Erasmus archives).

The main file is TI_Thesis_Draft.py. This includes the bulk of the code and requires grad.py, AggregateACS_DRAFT.py, and AdjacencyMatrix_DRAFT.py. 

Expendnm.prg, expendnm.run, auxpanel.dat, and panel.dat are required for GAUSS estimation. The code for specific models that can be plugged into expendnm.prg can be found in gmm_x.text and lfm_x.txt. To swap panel.dat to other .csv files with the same structure, use the code in load_in_routine.txt.

Neighbour_working.py is code to extract neighbour relationships from GIS. The code was tested in QGIS 3.40.4.

Given a working directory in use, the following subfolders should be made in order to make the code run smoothly. [x.csv] denote necessary datafiles that are not in GitHub, and [/x/] denotes zips of datafiles: <br/>
/WD/ <br/>
  TI_Thesis_Draft.py <br/>
  grad.py <br/>
  AggregateACS_DRAFT.py <br/>
  AdjacencyMatrix_DRAFT.py <br/>
  <br/>
  /expend/ <br/>
    Expendnm.prg  <br/>
    expendnm.run <br/>
    auxpanel.dat <br/>
    panel.dat <br/>
    <br/>
  /GIS/ <br/>
    neighbour_working.py <br/>
    [cb_2018_36_tract_500k.zip] <br/>
    [cb_2018_us_state_500k.zip] <br/>
    <br/>
  /Data/<br/>
  [inflation_data.csv]<br/>
  [Solar_Electric_Programs_Reported_by_NYSERDA__Beginning_2000_20241231.csv] <br/>
  /Census/<br/>
  [/B02001_NY_race/]<br/>
  [/DP05_NY_age/]<br/>
  [/S1101_NY_housing/]<br/>
  [/S1901_income/]<br/>
  [2020_2010_diff_cleaned.txt]<br/>
  <br/>
  /Simulations/<br/>
  <br/>
  /Plots/
  
  
  
  
  

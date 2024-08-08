# geostatistics_helper
Scripts for geostatistical analysis, EDA, kriging, ...


Design goal for optimising scripts was to reduce the total number of times the run button was pressed and collating/presenting results.

Calling the header file first ensures functions are preloaded and ready for use. 

Output maps are arranged for side-by-side comparison, with output data saved as csv in results folder.


Core functionality split across secondary scripts:
- A2_header.R   (loads libraries and other scripts)
- A2_functions.R   (Intermediate functions for EDA and different kriging methods)
- A2_optimised.R   (main script with user parameters, etc)


Current status is WIP (work in progresss), some issues still present and universal/indicator kriging functions not fully implemented.


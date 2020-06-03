# Synapse-Model
Short term synaptic plasticity models

## Dittman Models
These are models replicating the methods of Dittman et al. 2000:

Dittman JS, Kreitzer AC, Regehr WG. Interplay between facilitation, depression, and residual calcium at three presynaptic terminals. J Neurosci. 2000;20(4):1374‐1385. doi:10.1523/JNEUROSCI.20-04-01374.2000

**To run the model, simply download the desired model's .py file and call it using "python3 filename.py" in a bash interpreter. The necessary packages are found in the first few lines after the "import" keyword, these packages can be installed via pip/pip3 in your command line or from their respective source websites.**

*These files were written for python3, back compatibility with python2.7 or earlier is not certain*

___

## Two Pool and Maturation
This is a comparison of the model of synaptic vesicle release developed from the findings of Turcek et al. 2016 as well as my own maturation model:

Turecek et al., 2016, Cell Reports 17, 3256–3268
December 20, 2016
http://dx.doi.org/10.1016/j.celrep.2016.11.081

**To run the model, simply download Two_Pool_Model_func.m, Vesicle_maturation_model_func.m, and Synaptic_Model_Comparison_exported.m. Then open and run Synaptic_Model_Comparison_exported.m from Matlab.**

*Disclaimer: Two_Pool_Model_func.m is not the two pool model demonstrated in Turcek et al. 2016. It contains many of the same properties but appears to be missing some details. As such, the results using the parameters given in the Supplement of Turcek et al. will not match the results presented in their paper.*

The app should launch and bring you to a screen where you can freely change all parameters. You can switch between the maturation model and the two pool model with the tab at the top of the app. After entering the desired values for each parameter, click the "Set Parameters" button. The plot to the right and the table below should auto-update after pressing "Set Parameters".

If you want to plot new parameters, simply enter the new parameters and then click the "Set Parameters" button again. If you want to plot the variable which is currently selected in the dropdown, simply click the "Refresh" button and the plot should update. If you want to plot a different variable, simply select it from the dropdown menu and the plot should auto-update.

Some screenshots of how the app should look (on Windows) are available under [Synapse-Model/Examples](Examples)

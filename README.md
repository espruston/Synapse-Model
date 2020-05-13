# Synapse-Model
Synaptic facilitation models

This is a model of synaptic vesicle release developed from the findings of Turcek et al. 2016:

Turecek et al., 2016, Cell Reports 17, 3256–3268
December 20, 2016
http://dx.doi.org/10.1016/j.celrep.2016.11.081

**To run the model, simply download Two_Pool_Model_func.m, Vesicle_maturation_model_func.m, and Synaptic_Model_Comparison_exported.m. Then open and run Synaptic_Model_Comparison_exported.m from Matlab.**

*Disclaimer: Two_Pool_Model_func.m is not the two pool model demonstrated in Turcek et al. 2016. It contains many of the same properties but appears to be missing some details. As such, the results using the parameters given in the Supplement of Turcek et al. will not match the results presented in their paper.*

The app should launch and bring you to a screen where you can freely change all parameters. You can switch between the maturation model and the two pool model with the tab at the top of the app. After entering the desired values for each parameter, click the "Set Parameters" button. The plot to the right and the table below should auto-update after pressing "Set Parameters".

If you want to plot new parameters, simply enter the new parameters and then click the "Set Parameters" button again. If you want to plot the variable which is currently selected in the dropdown, simply click the "Refresh" button and the plot should update. If you want to plot a different variable, simply select it from the dropdown menu and the plot should auto-update.

Some screenshots of how the app should look (on Windows) are available under Synapse-Model/Examples

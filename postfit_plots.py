#!/usr/bin/env python
# coding: utf-8

"""
produce post fit plots using a pyhf fit in cabinetry

This uses a binned pyhf workspace as in dummy_plotting_spec.json

We fit in the two control regions, then extrapolate the fit into the
other region for a postfit plot

author : Daniel Noel

"""

import cabinetry
import numpy as np
import json
import pyhf
import os
import matplotlib.pyplot as plt



# these are how the spec has been defined
plot_info = {
           'pt1'   : {'PlotNBins' : 30, 'PlotBinLow' : 0., 'PlotBinHigh' : 200., 'titleX' : r'Leading $p_{T}$ [GeV]'},
           'pt2'   : {'PlotNBins' : 10, 'PlotBinLow' : 0., 'PlotBinHigh' : 200., 'titleX' : r'Subleading $p_{T}$ [GeV]'},
        }


def main():

    #First fit in the CRs only
    json_filepath="dummy_plotting_spec.json"
    with open(json_filepath) as serialized:
        spec = json.load(serialized)
    spec["measurements"][0]["config"]["poi"] = "" # needed for BkgOnly

    channel_list=["CR1", "CR2"]
    fit_results_main = fit_in_CRs(spec, channel_list)

    #get the full workspace
    with open(json_filepath) as serialized:
        spec = json.load(serialized)
    spec["measurements"][0]["config"]["poi"] = "" # needed for BkgOnly
    ws = pyhf.Workspace(spec)


    #list the variables to plot
    #this could be obtained from ws.channels, just picking the binned channels that correspond to variables
    vars_to_plot = ['pt1', 'pt2']
    region = 'CR'

    for variable in vars_to_plot:
        plot_variable(ws, fit_results_main, variable, vars_to_plot, region, plot_info)




def match_model(model, fit_results_main): 
    """
        Returns the model, postfit using the fit in CRs as stored in fit_results_main

        Matches results from a fit to a model by adding or removing parameters as needed.
        e.g. add dummy gamma factors for each bin

        Args
        ----------
        model : pyhf.pdf.Model
            Our initial model
        fit_results_main : cabinetry.fit.results_containers.FitResults
            Fit results in CRs

        Returns
        ----------
        model_postfit : pyhf.pdf.Model
            The postfit model
        """
    
  
    fit_results = cabinetry.model_utils.match_fit_results (model, fit_results_main)

    best_fit = fit_results[0]
    uncertainty = fit_results[1]
    corr_mat = fit_results[3]

    #look at post fit plots
    fit_results_matched = cabinetry.fit.results_containers.FitResults(
            np.asarray(best_fit), np.asarray(uncertainty), model.config.par_names(), corr_mat, 0
        )
    model_postfit = cabinetry.model_utils.prediction(model, fit_results=fit_results_matched)
    
    return model_postfit


def plot_datamc(model_postfit, data, region, var_name, plot_info):
    """
        Run the plotting using cabinetry.visualize.data_mc

        Args
        ----------
        model_postfit : pyhf.pdf.Model
            Our postfit model
        data: List
            The data in the region
        region : str
            Which region we're plotting
        var_name : str
            Which variable we're plotting
        plot_info: dict
            Dictionary of plotting information

        """
    #plot the postfit result
    cabinetry.visualize.data_mc(model_postfit, data,log_scale=True)
    print("plotting, ", plot_info[var_name])
    
    fig = plt.gcf()
    ax1, ax2 = fig.axes

    ax2.set_ylabel("data / SM")
    ax1.set_ylabel("Events")
    ax1.yaxis.set_label_coords(-0.12,0.95) # move "Events"
    plt.xlabel(plot_info[var_name]['titleX'])

    # change ticklabels to our data
    # note the data has bins from 0->nbins
    nticks = 6
    ticklabels = np.linspace(plot_info[var_name]['PlotBinLow'], plot_info[var_name]['PlotBinHigh'], nticks)
    ticklabels = [float('{:.2f}'.format(i)) for i in ticklabels]
    tick_locs = np.linspace(min(ax2.get_xticks()), max(ax2.get_xticks()), nticks )
    ax2.set_xticks(tick_locs)
    ax2.set_xticklabels(ticklabels)

    
    for txt in fig.texts: txt.set_visible(False) # remove text on there currently
    
    out_dir = os.path.join('plots',region)
    os.makedirs(out_dir, exist_ok=True)
    plt.savefig(os.path.join(out_dir,"{}_postfit.png".format(var_name)), transparent=False)
    plt.savefig(os.path.join(out_dir,"{}_postfit.pdf".format(var_name)), transparent=False)

    plt.clf()
    plt.close('all')


def plot_variable(ws, fit_results_main, variable, all_variables, region, plot_info):
    """
        Run the plotting on a variable, first getting the postfit yields,
        then running the plotting function

        Args
        ----------
        model_postfit : pyhf.pdf.Model
            Our postfit model
        fit_results_main : cabinetry.fit.results_containers.FitResults
            Fit results in CRs
        variable : str
            Which variable we're plotting
        all_variables : List
            List of all the variables in the region
        region : str
            Which region we're plotting
        plot_info: dict
            Dictionary of plotting information

        """
    print(vars(ws))
    #remove variables we are not plotting
    variables_to_prune = all_variables[:]
    variables_to_prune.remove(variable)
    ws_pruned = pyhf.Workspace.prune(ws,channels= variables_to_prune)
    model, data = cabinetry.model_utils.model_and_data(ws_pruned)
    
    print(variable, variables_to_prune)

    #get postfit model
    model_postfit = match_model(model, fit_results_main)
    
    #the way it's saved we can extract the variable like this:
    var_name = variable.replace(region+"_", "")
    
    #plot variable "var_name" in region "region"
    plot_datamc(model_postfit, data, region, var_name, plot_info)


def fit_in_CRs(spec, CR_list):
    """ Fit just in the CRs

        Args
        ----------
        spec : dict
            json spec of theinitial model

        Returns
        ----------

        fit_results_main : cabinetry.fit.results_containers.FitResults
            Fit results in CRs
            """
    spec["measurements"][0]["config"]["poi"] = "" # needed for BkgOnly

    # fit in the CRs only
    spec["channels"] = [c for c in spec["channels"] if c["name"] in CR_list]

    ws = pyhf.Workspace(spec)

    #fit the model in the CRs
    model, data = cabinetry.model_utils.model_and_data(ws)
    fit_results_main = cabinetry.fit.fit(model, data)

    return fit_results_main


if __name__ == "__main__":
    main()


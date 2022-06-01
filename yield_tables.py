#!/usr/bin/env python
# coding: utf-8
"""
Writing yield tables in Cabinetry, using a fit in CRs and
extrapolating them to VRs and SRs


author : Daniel Noel
"""
import json
import os
import pickle

import cabinetry
import numpy as np
import pyhf
from tabulate import tabulate


def main():

    # first fit the model in the CRs
    json_filepath = "dummy_spec.json"

    with open(json_filepath) as serialized:
        spec = json.load(serialized)

    spec["measurements"][0]["config"]["poi"] = ""  # needed for BkgOnly
    channel_list = ["CR1", "CR2"]  # only fit CRs
    spec["channels"] = [c for c in spec["channels"] if c["name"] in channel_list]
    ws = pyhf.Workspace(spec)

    # fit the model
    model, data = cabinetry.model_utils.model_and_data(ws)
    fit_results_main = cabinetry.fit.fit(model, data)

    # ----------------
    # Making the table
    # reopen the model and apply the previous fit to the other regions to generate the yield table
    with open(json_filepath) as serialized:
        spec = json.load(serialized)
    spec["measurements"][0]["config"]["poi"] = ""  # needed for BkgOnly

    # and manually edit the tex (todo edit them in here e.g. add the hlines)

    ws = pyhf.Workspace(spec)
    yt = yield_table(ws, fit_results_main)
    yt.write_yield_table(4)
    yt.save_histfitter_table()


def get_channels_to_prune(spec, sample):
    """pyhf.Workspace.prune() won't work to leave an totally empty channel
    so if pruning away to one background, must only include channels where the background exists
    i.e. prune off empty channels too
    
    
    Args
    ----------
    spec : dict
        The json spec for input to pyhf
    sample : str
        The samcorr_matple we wish to use
        
    Returns
    ----------
    channels_to_prune : List
        List of the channels where the sample is not present,
        so need to be pruned off when everything but the sample
        is being pruned
    """
    channels_to_prune = []
    for channel in spec["channels"]:
        # sample that exist in the channel
        sample_names = [
            channel["samples"][i]["name"] for i in range(len(channel["samples"]))
        ]

        # if our samplele isn't in the channel we need to remove it!
        if sample not in sample_names:
            channels_to_prune.append(channel["name"])

    return channels_to_prune


def match_fit_and_get_yields_stdev(ws, fit_results):
    """Using fit_results from a previous fit get the post-fit 
    yield and errors for a given ws. 
    Note - finds overall combined stdev of all samples in ws. 
    So this ws can be pruned accordingly before being passed into 
    this function in order to get individual sample yields.
    
    
    Args
    ----------
    ws : pyhf.Workspace
        The workspace we wish to find the yields from
    fit_results : cabinetry.fit.results_containers.FitResults
        The fitted parameters from a previous fit:
        fit_results = cabinetry.fit.fit(model, data)
        
    Returns
    ----------
    yields : np.array
        List of yields for each channel in ws
    stdev : List
        List of standard deviations for each channel in ws
        
    """

    # match the fit results to the previous fit using cabinetry
    model, data = cabinetry.model_utils.model_and_data(ws)
    fit_results_matched = cabinetry.model_utils.match_fit_results(model, fit_results)
    best_fit = fit_results_matched[0]
    uncertainty = fit_results_matched[1]
    corr_mat = fit_results_matched[3]

    # use the already fit best_fit, uncertainty,corr_mat
    stdev = cabinetry.model_utils.yield_stdev(model, best_fit, uncertainty, corr_mat)

    # verify the std output are the same (no binned channels)
    for i in range(len(stdev[1])):
        assert stdev[0][i][0] == stdev[1][i]

    yields = model.expected_data(best_fit, include_auxdata=False)
    return yields, stdev[1]


def get_sample_yields_stdev(ws, fit_results_prev):
    """Using fit_results from a previous fit get the post-fit yields and 
    errors for each individual sample in each different region
    
    Typical use case: Fit in the CRs, use the fit_results from CRs
    to calculate post-fit yields and errors in VRs and SRs
    
    Note that due to correlations, the errors must be calculated for each sample
    individually and may not necessarily sum in quadrature to the overall error
    on the combined yield from all the samples
    
    Args
    ----------
    ws : pyhf.Workspace
        The workspace we wish to find the yields from
    fit_results_prev : cabinetry.fit.results_containers.FitResults
        The fitted parameters from a previous fit:
        fit_results = cabinetry.fit.fit(model, data)
        
    Returns
    ----------
    postfit_sample_yields : dict
        dictionary of yields for each sample
        {sample : {region : yield }}
    postfit_sample_stdev : dict
        dictionary of stdev for each sample
        {sample : {region : stdev }}
    """
    model = ws.model()  # Note can be quicker to run once before
    spec = model.spec

    # initialise dictionaries
    postfit_sample_stdev = {}
    postfit_sample_yields = {}

    # calculate yield and stdev for each sample individually
    all_samples = ws.samples
    for sample in all_samples:

        # in order to calculate yield and stdev for a single sample
        # we must prune off all other samples
        samples_to_prune = all_samples[:]
        samples_to_prune.remove(sample)
        print(sample, samples_to_prune)

        # channels where the sample is not present
        channel_list = get_channels_to_prune(spec, sample)

        # remove all samples apart from sample
        ws_pruned = pyhf.Workspace.prune(
            ws, channels=channel_list, samples=samples_to_prune
        )
        print(ws_pruned.samples)

        yields, stdev = match_fit_and_get_yields_stdev(ws_pruned, fit_results_prev)

        # save the yields and stdev
        postfit_sample_stdev[sample] = dict(zip(ws_pruned.channels, stdev))
        postfit_sample_yields[sample] = dict(zip(ws_pruned.channels, yields))

    return postfit_sample_yields, postfit_sample_stdev


def get_prefit_fit_results(model):
    """cabinetry.model_utils.prefit_uncertainties gives 0 for all
    https://github.com/scikit-hep/cabinetry/blob/a79410a4fcdb97d0cc797cce8806f6d9780c84be/src/cabinetry/model_utils.py#L125
    
    It seems we need to give different uncertainties for mu and lumi parameters.
    Let's rebuild everything from scratch...
    Also this function packages up into a fit_results container ready for matching 
    to a workspace
    
    
    Args
    ----------
    model : pyhf.pdf.Model
        The model we wish to find the prefit parameters from
        
    Returns
    ----------
    fit_results : cabinetry.fit.results_containers.FitResults
        The prefit fit parameters
        fit_results = cabinetry.fit.fit(model, data)
    """

    uncertainty = []
    best_fit = []
    for i in range(model.config.npars):

        # stat error and mu params
        if (
            "staterror" in model.config.parameters[i]
            or "mu" in model.config.parameters[i]
        ):
            best_fit.append(1.0)
            uncertainty.append(0.0)
        elif "lumi" in model.config.parameters[i]:
            best_fit.append(1.0)
            # this can be found from model.config.suggested_bounds()
            # but for now hard code lumi uncertainty
            uncertainty.append(0.039)
        else:
            # for the systematics
            best_fit.append(0.0)
            uncertainty.append(1.0)

    uncertainty = np.array(uncertainty)
    best_fit = np.array(best_fit)
    corr = np.identity(model.config.npars)  # identity matrix for correlation prefit

    from cabinetry.fit.results_containers import FitResults

    fit_results_prefit = FitResults(
        np.asarray(best_fit),
        np.asarray(uncertainty),
        model.config.parameters,
        corr,
        0,  # best_twice_nll
        0,  # goodness_of_fit
    )

    return fit_results_prefit


def round_yields(*args):
    """Takes in yields and rounds them according to the pdg recommendations:
    https://pdg.lbl.gov/2011/reviews/rpp2011-rev-rpp-intro.pdf page 13
    
    
    Args
    ----------
    *args : float
        numbers to be rounded
    
    """
    smallest_yield = float(min(args))

    # The recommended rounding convention
    if smallest_yield > 3.55:
        number_dp = 0
    elif smallest_yield > 0.355:
        number_dp = 1
    elif smallest_yield > 0.0355:
        number_dp = 2
    elif smallest_yield > 0.00355:
        number_dp = 3
    elif smallest_yield > 0.000355:
        number_dp = 4
    elif smallest_yield > 0.0000355:
        number_dp = 5
    elif smallest_yield == 0:
        number_dp = 2
    else:
        raise Exception("Issue - none of the yields should be negative")

    # round(arg,0) gave .0 at the end, we don't want this
    if number_dp == 0:
        rounded_yields_list = [round(arg) for arg in args]
    else:
        rounded_yields_list = [round(arg, number_dp) for arg in args]

    # return just a number in case of one argument
    if len(args) == 1:
        return rounded_yields_list[0]
    else:
        return rounded_yields_list


# def round_yields(*args):
#     """If want to use 2dp throughout
#     """
#     #2dp for all
#     rounded_yields_list = [round(arg,2) for arg in args]

#     #return just a number in case of one argument
#     if len(args) == 1: return rounded_yields_list[0]
#     else: return rounded_yields_list


def get_yield_strings(yields, stdev):
    """
    Get the strings for writing in the table:
    from the mean and up/down variation of yields obtain a latex string of the yield
    
    
    Args
    ----------
    yield : Dict
        The total yields per channel
        {channel : yield}
    stdev: Dict
        The total errors per channel
        {channel : stdev}
    Returns
    ----------
    strings : Dict
        contains strings for writing in the table
        in latex format 
       {channel : string }
    """

    strings = {}

    # for each channel get the latex yield string
    for channel in yields:
        mean_sample_yield = float(yields[channel])
        std = stdev[channel]
        print(mean_sample_yield, std)

        # TODOO check this may affect the error
        if mean_sample_yield < 0:
            mean_sample_yield = 0.01  # round up to 0.01 if negative for FNP
        if mean_sample_yield < 0.01:
            mean_sample_yield = 0.0  # round other bkgs to 0

        up_variation = std
        if mean_sample_yield - std < 0:
            down_variation = mean_sample_yield  # truncate at 0
        else:
            down_variation = std

        print(channel, mean_sample_yield, up_variation, down_variation)
        # rounded other bkgs to 0
        if mean_sample_yield == 0.0:
            up_variation = 0.0
            down_variation = 0.0

        (mean_sample_yield, up_variation, down_variation) = round_yields(
            mean_sample_yield, up_variation, down_variation
        )
        # symmetric
        if up_variation == down_variation:
            strings[channel] = "${0} \pm {1}$".format(mean_sample_yield, up_variation)
        # asymmetric (around 0)
        else:
            strings[channel] = "${0}^{{+{1}}}_{{-{2}}}$".format(
                mean_sample_yield, up_variation, down_variation
            )
    return strings


def get_sample_yield_strings(sample_yields, sample_stdev):
    """Get the strings for writing in the table:
    from the mean and up/down variation of yields obtain a latex string of the yield
    
    
    Args
    ----------
    sample_yields : Dict
        The sample yields per channel
        {sample : {channel : yield} }
    sample_stdev: Dict
        The sample errors per channel
        {sample : {channel : stdev}}
    Returns
    ----------
    sample_strings : Dict
        strings for writing in the table
        {sample : {channel : string}}
        
    """
    sample_strings = {}
    for sample in sample_yields:

        sample_strings[sample] = {}
        yields = sample_yields[sample]
        stdev = sample_stdev[sample]

        # get yield strings for each sample
        sample_strings[sample] = get_yield_strings(yields, stdev)

    return sample_strings


class yield_table:
    """A class to collect yields and write latex yield tables
    
    Example usage : write a yield table with 4 columns
    ----------
        yt = yield_table(ws,fit_results_main)
        yt.write_yield_table(4)

    """

    def __init__(self, ws, fit_results_main):
        """
        Collect the relevant yields
        
        Args
        ----------
        ws : pyhf.Workspace
            The workspace we wish to find the yields from

        fit_results_main : cabinetry.fit.results_containers.FitResults
            The fitted parameters from a previous fit:
            fit_results = cabinetry.fit.fit(model, data)

        Usage
        """
        self.ws = ws
        self.fit_results_main = fit_results_main

        self.all_samples = ws.samples
        self.channels = ws.channels
        self.row_names = self.get_row_names()

        model, data = cabinetry.model_utils.model_and_data(ws, include_auxdata=False)
        self.observed_data = dict(zip(self.channels, data))

        # _______________________________________
        # prefit
        self.fit_results_prefit = get_prefit_fit_results(model)

        # prefit sample yields
        self.prefit_sample_yields, self.prefit_sample_stdev = get_sample_yields_stdev(
            self.ws, self.fit_results_prefit
        )

        # prefit overall
        self.prefit_yields, self.prefit_stdev = match_fit_and_get_yields_stdev(
            ws, self.fit_results_prefit
        )
        self.prefit_stdev = dict(zip(self.channels, self.prefit_stdev))
        self.prefit_yields = dict(zip(self.channels, self.prefit_yields))

        # _______________________________________
        # postfit sample yields
        self.postfit_sample_yields, self.postfit_sample_stdev = get_sample_yields_stdev(
            self.ws, self.fit_results_main
        )

        # postfit overall
        self.postfit_yields, self.postfit_stdev = match_fit_and_get_yields_stdev(
            self.ws, self.fit_results_main
        )
        self.postfit_stdev = dict(zip(ws.channels, self.postfit_stdev))
        self.postfit_yields = dict(zip(ws.channels, self.postfit_yields))

        # _______________________________________
        # def get strings for table
        self.prefit_strings = get_yield_strings(self.prefit_yields, self.prefit_stdev)
        self.prefit_samples_strings = get_sample_yield_strings(
            self.prefit_sample_yields, self.prefit_sample_stdev
        )

        self.postfit_strings = get_yield_strings(
            self.postfit_yields, self.postfit_stdev
        )
        self.postfit_samples_strings = get_sample_yield_strings(
            self.postfit_sample_yields, self.postfit_sample_stdev
        )

    def get_row_names(self):
        """get a list of the row names:
           Observed data
           Fitted yield (total)
           individual fitted background yields
           prefit yield (total)
           individual prefit background yields"""
        row_names = (
            ["Observed data", "Fitted yield"]
            + ["Fitted {} events".format(samp) for samp in self.all_samples]
            + ["Prefit yield"]
            + ["MC exp. {} events".format(samp) for samp in self.all_samples]
        )
        return row_names

    def get_table_content_dict(self, NCOLS):
        """
        Get the dictinary for writing the table with tabulate
        
        Args
        ----------
        NCOLS : int
            The number of columns per page to plot
        
        Returns
        ----------
        table_dict : Dict
            dictionary for writing table with tabulate
            
            {channel name : List } 
            List of the column content 
        """
        table_dict = {}
        for i, channel_name in enumerate(self.channels):

            channel_name_latex = channel_name.replace("_", " ").replace(" cuts", "")
            table_dict[channel_name_latex] = []

            # add observed data
            table_dict[channel_name_latex].append(int(self.observed_data[channel_name]))

            # add fitted total yields
            table_dict[channel_name_latex].append(self.postfit_strings[channel_name])

            # add postfit yields for each sample
            for sample in self.all_samples:
                if channel_name in self.postfit_samples_strings[sample]:
                    table_dict[channel_name_latex].append(
                        self.postfit_samples_strings[sample][channel_name]
                    )
                # sample not in the region
                else:
                    table_dict[channel_name_latex].append("$0.00 \pm 0.00$")

            # add total prefit yields
            table_dict[channel_name_latex].append(self.prefit_strings[channel_name])

            # add prefit yields for each sample
            for sample in self.all_samples:
                if channel_name in self.postfit_samples_strings[sample]:
                    table_dict[channel_name_latex].append(
                        self.prefit_samples_strings[sample][channel_name]
                    )
                # sample not in the region
                else:
                    table_dict[channel_name_latex].append("$0.00 \pm 0.00$")

            # return only NCOLS at a time
            if i % NCOLS + 1 == NCOLS:
                yield table_dict
                table_dict = {}

        # flush the final table dict out
        if len(table_dict) != 0:
            yield table_dict

    def write_yield_table(self, NCOLS, name=""):
        """Write the latex yield table to file, and call pdflatex to
        generate the pdf

        Args
        ----------
        NCOLS : int
            The number of columns per page to plot
        """

        print("Writing Table...")
        file_name = "YieldTable{}.tex".format(name)
        f = open(file_name, "w")
        f.write("\\documentclass{article}\n")
        f.write("\\usepackage[a4paper,margin=1in,landscape]{geometry}\n")
        f.write("\\usepackage{booktabs}\n")
        f.write("\\usepackage{graphicx}\n")
        f.write("\\begin{document}\n")

        for table_dict in self.get_table_content_dict(NCOLS):
            print(table_dict, table_dict.keys(), self.row_names)
            f.write(
                tabulate(
                    table_dict,
                    headers=table_dict.keys(),
                    showindex=self.row_names,
                    tablefmt="latex_raw",
                )
            )
            f.write("\n\n")

        f.write("\\end{document}\n")

        f.close()

        os.makedirs("output", exist_ok=True)
        os.rename(file_name, "output/{0}".format(file_name))

    def save_histfitter_table(self):
        """get a pickle in the histfitter format for use with
        PullPlot.py scripts"""
        histfitter_table = {
            "names": self.channels,
            "nobs": list(self.observed_data.values()),
            "TOTAL_MC_EXP_BKG_events": list(self.prefit_yields.values()),
            "TOTAL_MC_EXP_BKG_err": list(self.prefit_stdev.values()),
            "TOTAL_FITTED_bkg_events": list(self.postfit_yields.values()),
            "TOTAL_FITTED_bkg_events_err": list(self.prefit_stdev.values()),
        }

        # sample specific
        for name, yield_dict in {
            "Fitted_events": self.postfit_sample_yields,
            "Fitted_err": self.postfit_sample_stdev,
            "MC_exp_events": self.prefit_sample_yields,
            "MC_exp_err": self.prefit_sample_stdev,
        }.items():
            for bkg in yield_dict:
                postfit_yields = []
                for ch in self.channels:
                    if ch in yield_dict[bkg]:
                        postfit_yields.append(yield_dict[bkg][ch])
                    else:
                        postfit_yields.append(0.0)
                histfitter_table["{}_{}".format(name, bkg)] = postfit_yields

        os.makedirs("output", exist_ok=True)
        with open("output/yield_table_pyhf.pickle", "wb") as handle:
            pickle.dump(
                histfitter_table, handle, protocol=2
            )  # if we want to use python2 for the HF use this protocol...


if __name__ == "__main__":
    main()

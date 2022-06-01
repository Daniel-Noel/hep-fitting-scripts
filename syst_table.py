#!/usr/bin/env python
# coding: utf-8

"""
Output systematic breakdowns to json
as method 1. in histfitter paper https://arxiv.org/pdf/1410.1280.pdf#page=27

author : Daniel Noel
"""
import json
import os
import time

import cabinetry
import pyhf
from cabinetry.fit.results_containers import FitResults

from yield_tables import match_fit_and_get_yields_stdev


def main():
    json_filepath = "dummy_spec.json"
    i_ch = 2  # index of channel interested in

    # do a fit in the CR without stat errors to get total syst error component
    ws = get_ws_CR_only(json_filepath, no_stat=True)
    model, data = cabinetry.model_utils.model_and_data(ws)
    fit_results_main = cabinetry.fit.fit(model, data)

    # calculate the total uncert
    ws = get_ws_allchannels(json_filepath)
    yields_full, stdev_full = match_fit_and_get_yields_stdev(ws, fit_results_main)
    stdevs = {"full": stdev_full[i_ch]}

    # calculate contribution of each systematic one by one
    labels = model.config.parameters
    print(labels)
    for syst in labels:
        _, stdev = get_stdev_syst(ws, fit_results_main, [syst])
        stdevs[syst] = stdev[i_ch]

    # consider systematics grouped together
    grouped_systs = ["syst_1", "syst_2"]
    _, stdev = get_stdev_syst(ws, fit_results_main, grouped_systs)
    stdevs["combined"] = stdev[i_ch]

    # calculate the percentage uncertainty too
    perc_uncerts = {}
    for syst in stdevs:
        perc_uncerts[syst] = stdevs[syst] / yields_full[i_ch] * 100

    # save the results
    os.makedirs("output", exist_ok=True)
    with open(os.path.join("output", "syst_uncerts.json"), "w+") as out_file:
        json.dump(stdevs, out_file, indent=4)
    with open(os.path.join("output", "syst_percent_uncerts.json"), "w+") as out_file:
        json.dump(perc_uncerts, out_file, indent=4)


def get_ws_CR_only(json_filepath, no_stat=False):
    """no stat uncerts"""
    with open(json_filepath) as serialized:
        spec = json.load(serialized)

    spec["measurements"][0]["config"]["poi"] = ""  # needed for BkgOnly
    channel_list = ["CR1", "CR2"]  # only fit CRs
    spec["channels"] = [c for c in spec["channels"] if c["name"] in channel_list]
    ws = pyhf.Workspace(spec)

    # todo - check if works and remove this / make option
    if no_stat:
        return ws.prune(["staterror_CR1", "staterror_CR2"])
    else:
        return ws


def get_ws_allchannels(json_filepath):
    with open(json_filepath) as serialized:
        spec = json.load(serialized)
    spec["measurements"][0]["config"]["poi"] = ""  # needed for BkgOnly
    ws = pyhf.Workspace(spec)
    return ws


def get_stdev_syst(ws, fit_results, systs):
    """using a full workspace ws, and previous fit fit_results,
    we calculate the stdev due to the systematics in the list systs """

    # get a workspace with just desired systematics for the error calculation
    labels = ws.model().config.parameters
    for syst in systs:
        assert syst in labels  # if not then provide correct named systs
        labels.remove(syst)
    ws = ws.prune(labels)

    yields, stdevs = match_fit_and_get_yields_stdev(ws, fit_results)

    return yields, stdevs


if __name__ == "__main__":
    main()

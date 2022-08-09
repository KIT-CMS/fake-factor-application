#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import yaml
import ROOT
import numpy
import copy
from array import array
import six
from numpy import sqrt, power
import logging
logger = logging.getLogger()


def determine_fractions(config, configkey, fractionsfile, era, emb=True, tag='smhtt', era_labels=None):
    configpath = copy.deepcopy(config)
    if isinstance(configpath, six.string_types) and len(configpath):
        config = yaml.load(open(configpath))
        if "categories" not in config.keys():
            logger.critical("Config file %s has to contain key 'categories'!" % configpath)
            raise Exception
    else:
        assert isinstance(configpath, dict), "Neither config file nor config list provided"

    categories = config["categories"]
    hist_file = ROOT.TFile(fractionsfile, "READ")
    if era_labels is None:
        era_labels = {
            "2016": "Run2016",
            "2017": "Run2017ReReco31Mar"
        }
    composition = {
        "mt": {
            "data": ["data_obs"],
            "W": ["W", "VVJ", "ZJ"],
            "TT": ["TTJ"],
            "QCD": ["QCD"],
            "real": ["ZL", "TTL", "VVL"] + (["EMB"] if emb else ["ZTT", "TTT", "VVT"])
        },
        "et": {
            "data": ["data_obs"],
            "W": ["W", "VVJ", "ZJ"],
            "TT": ["TTJ"],
            "QCD": ["QCD"],
            "real": ["ZL", "TTL", "VVL"] + (["EMB"] if emb else ["ZTT", "TTT", "VVT"])
        },
        "tt": {
            "data": ["data_obs"],
            "W": ["W", "VVJ", "ZJ"],
            "TT": ["TTJ"],
            "QCD": ["QCD"],
            "real": ["ZL", "TTL", "VVL"] + (["EMB"] if emb else ["ZTT", "TTT", "VVT"])
        }
    }

    fractions = {}
    for channel in categories.keys():
        subdict = {}
        for category in categories[channel]:
            subsubdict = {}
            for fraction in composition[channel].keys():
                print("Determining fractions: #{ch}#{ch}_{cat}#data_obs#{tag}#{era}#{expr}#125#".format(
                    ch=channel,
                    cat=category,
                    tag=tag,
                    era=era_labels[era],
                    expr=configkey)
                )
                subsubdict[fraction] = copy.deepcopy(
                    hist_file.Get(
                        "#{ch}#{ch}_{cat}#data_obs#{tag}#{era}#{expr}#125#".
                        format(
                            ch=channel,
                            cat=category,
                            tag=tag,
                            era=era_labels[era],
                            expr=configkey)))
                subsubdict[fraction].Reset()
            for fraction in composition[channel].keys():
                if fraction == "QCD":
                    continue
                for process in composition[channel][fraction]:
                    subsubdict[fraction].Add(
                        hist_file.Get(
                            "#{ch}#{ch}_{cat}#{proc}#{tag}#{era}#{expr}#125#".
                            format(
                                ch=channel,
                                cat=category,
                                proc=process,
                                tag=tag,
                                era=era_labels[era],
                                expr=configkey)))
                if fraction == "data":
                    subsubdict["QCD"].Add(subsubdict[fraction], 1.0)
                else:
                    subsubdict["QCD"].Add(subsubdict[fraction], -1.0)
            # Normalize to data to get fractions
            denominator_hist = copy.deepcopy(subsubdict["data"])
            for fraction in composition[channel].keys():
                subsubdict[fraction].Divide(denominator_hist)
            # Sanity check: if QCD negative i.e. data < MC, normalize to MC
            for i in range(subsubdict["QCD"].GetNbinsX() + 2):
                qcd_fraction = subsubdict["QCD"].GetBinContent(i)
                if qcd_fraction < 0.0:
                    logger.info("Found bin with negative QCD fraction (%s, %s, index %i). Set fraction to zero and rescale other fractions..." % (channel, category, i))
                    subsubdict["QCD"].SetBinContent(i, 0)
                    for fraction in composition[channel].keys():
                        if not fraction == "data":
                            subsubdict[fraction].SetBinContent(i, subsubdict[fraction].GetBinContent(i) / (1.0 - qcd_fraction))
                            print("Rescaled %s fraction to %f" % (fraction, subsubdict[fraction].GetBinContent(i)))
            subdict[category] = subsubdict
        fractions[channel] = subdict
    hist_file.Close()
    return fractions


def apply_fake_factors(
        datafile, friendfilelists, outputfile, category_mode, fakefactordirectories,
        fractions, use_fractions_from_workspace, configpath, expression, era,
        pipeline_selection=None, treename="ntuple", eventrange=None, rootfilemode="update",
        use_mva_tauid=False,
):

    config = yaml.load(open(configpath))
    if "categories" not in config.keys():
        logger.critical("Config file %s has to contain key 'categories'!" % configpath)
        raise Exception
    categories = config["categories"]

    rootfilemode = rootfilemode.lower()
    if rootfilemode not in ["update", "recreate"]:
        logger.critical("Mode %s not appropriate for create of ROOT file. Please choose from 'update' and 'recreate'"%rootfilemode)
        raise Exception

    if pipeline_selection is not None:
        assert isinstance(pipeline_selection, six.string_types), "Invalid pipeline_selection"

    # documented in https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauJet2TauFakes
    unc_shifts = {
        "et": [
            "ff_qcd_syst", "ff_qcd_dr0_njet0_stat", "ff_qcd_dr0_njet1_stat", "ff_qcd_dr0_njet2_stat",
            "ff_w_syst", "ff_w_dr0_njet0_stat", "ff_w_dr0_njet1_stat", "ff_w_dr0_njet2_stat",
            "ff_w_dr1_njet0_stat", "ff_w_dr1_njet1_stat", "ff_w_dr1_njet2_stat",
            "ff_tt_syst", "ff_tt_stat", "ff_tt_morphed", "ff_tt_sf", "ff_corr_tt_syst",
            "ff_frac_w",
            "ff_w_lepPt",
            "ff_w_mc",
            "ff_w_mt",
            "ff_corr_w_lepPt",
            "ff_corr_w_mt",
            "ff_qcd_mvis",
            "ff_qcd_muiso",
            "ff_corr_qcd_mvis_osss",
            "ff_qcd_mvis_osss",
            "ff_corr_qcd_mvis",
            "ff_corr_qcd_muiso",
            "ff_qcd_mc",
            "ff_qcd_dr0_njet0_morphed_stat", "ff_qcd_dr0_njet1_morphed_stat", "ff_qcd_dr0_njet2_morphed_stat",
            "ff_w_dr0_njet0_morphed_stat", "ff_w_dr0_njet1_morphed_stat", "ff_w_dr0_njet2_morphed_stat",
            "ff_w_dr1_njet0_morphed_stat", "ff_w_dr1_njet1_morphed_stat", "ff_w_dr1_njet2_morphed_stat",
            "ff_tt_dr0_njet0_morphed_stat", "ff_tt_dr0_njet1_morphed_stat",
        ],
        "mt": [
            "ff_qcd_syst", "ff_qcd_dr0_njet0_stat", "ff_qcd_dr0_njet1_stat", "ff_qcd_dr0_njet2_stat",
            "ff_w_syst", "ff_w_dr0_njet0_stat", "ff_w_dr0_njet1_stat", "ff_w_dr0_njet2_stat",
            "ff_w_dr1_njet0_stat", "ff_w_dr1_njet1_stat", "ff_w_dr1_njet2_stat",
            "ff_tt_syst", "ff_tt_stat", "ff_tt_morphed", "ff_tt_sf", "ff_corr_tt_syst",
            "ff_frac_w",
            "ff_w_lepPt",
            "ff_w_mc",
            "ff_w_mt",
            "ff_corr_w_lepPt",
            "ff_corr_w_mt",
            "ff_qcd_mvis",
            "ff_qcd_muiso",
            "ff_corr_qcd_mvis_osss",
            "ff_qcd_mvis_osss",
            "ff_corr_qcd_mvis",
            "ff_corr_qcd_muiso",
            "ff_qcd_mc",
            "ff_qcd_dr0_njet0_morphed_stat", "ff_qcd_dr0_njet1_morphed_stat", "ff_qcd_dr0_njet2_morphed_stat",
            "ff_w_dr0_njet0_morphed_stat", "ff_w_dr0_njet1_morphed_stat", "ff_w_dr0_njet2_morphed_stat",
            "ff_w_dr1_njet0_morphed_stat", "ff_w_dr1_njet1_morphed_stat", "ff_w_dr1_njet2_morphed_stat",
            "ff_tt_dr0_njet0_morphed_stat", "ff_tt_dr0_njet1_morphed_stat", 
        ],
        "tt": [
            "ff_qcd_syst", "ff_qcd_dr0_njet0_stat", "ff_qcd_dr0_njet1_stat", "ff_qcd_dr0_njet2_stat",
            "ff_w_syst", "ff_tt_syst", "ff_w_frac_syst", "ff_tt_frac_syst",
            "ff_qcd_mvis",
            "ff_corr_qcd_mvis_osss",
            "ff_qcd_mvis_osss",
            "ff_qcd_tau2_pt",
            "ff_corr_qcd_mvis",
            "ff_corr_qcd_tau2_pt",
            "ff_qcd_mc",
            "ff_qcd_dr0_njet0_morphed_stat", "ff_qcd_dr0_njet1_morphed_stat", "ff_qcd_dr0_njet2_morphed_stat",
        ]
    }

    # Determine channel
    channels = ["et", "mt", "tt"]

    selected_channels = [channel for channel in fractions.keys() if channel in channels]
    if len(selected_channels) > 1:
        raise Exception("More than one channel configured, which is not forseen !")
    elif len(selected_channels) == 0:
        print('done (no matching channels found)')
        return 1

    channel = selected_channels[0]
    input_file = ROOT.TFile(datafile, "READ")
    input_tree = input_file.Get(treename)
        # Load fake factors histograms
    ff_file = ROOT.TFile.Open(fakefactordirectories[channel])
    print(fakefactordirectories[channel])
    ff = ff_file.Get('ff_comb')

    # Prepare output
    print("...initialize %s" % outputfile)
    if not os.path.exists(os.path.dirname(outputfile)):
        os.mkdir(os.path.dirname(outputfile))
    output_file = ROOT.TFile(outputfile, rootfilemode)
    # output_root_dir = output_file.mkdir("%s_%s" % (channel, pipeline))
    # output_root_dir.cd()
    output_tree = ROOT.TTree(treename, treename)
    # determine all variables, and all the corresponding shifts available, in order to run all those shifts as piplines
    default_variables = []
    if channel == "et":
        default_variables = ["pt_1", "pt_2", "deltaR_ditaupair", "m_vis", "njets", "iso_1", "mt_1_pf", "mt_1", "decaymode_2"]
    if channel == "mt":
        default_variables = ["pt_1", "pt_2", "deltaR_ditaupair", "m_vis", "njets", "iso_1", "mt_1_pf", "mt_1", "decaymode_2"]
    elif channel == "tt":
        default_variables = ["pt_1", "pt_2", "deltaR_ditaupair", "m_vis", "njets", "iso_1", "mt_1_pf", "mt_1", "decaymode_1", "decaymode_2"]
    vardict= {}

    quantites = [x.GetName() for x in input_tree.GetListOfLeaves()]
    for variable in default_variables:
        if variable not in vardict.keys():
            vardict[variable] = set()
        for quant in quantites:
            if quant.split("__")[0] == variable:
                if len(quant.split("__")) > 1:
                   vardict[variable].add(quant.split("__")[1])
    required_piplines = list(set([x for var in vardict.keys() for x in vardict[var]])) + ["nominal"]
    # remove pipelines with jes in name from the list of required pipelines
    required_piplines = [x for x in required_piplines if "jes" not in x and "jer" not in x and "met" not in x]
    print("Running pipelines: {}".format(required_piplines))
    output_buffer = {}
    varlist_dict = {}

    # for channel in selected_channels:
    print("Channel: {}".format(channel))
    for pipeline in required_piplines:

        # if it is nominal, just use the default variable, if it is a shift, try to find, if the shifted variable exists, use it, otherwise use the default
        if pipeline == "nominal":
            varlist_dict["nominal"] = []
            for var in default_variables:
                if var in quantites:
                    varlist_dict["nominal"].append(var)
        else:
            varlist_dict[pipeline] = []
            for var in default_variables:
                if "{var}__{pipeline}".format(var=var, pipeline=pipeline) in quantites:
                    varlist_dict[pipeline].append("{var}__{pipeline}".format(var=var, pipeline=pipeline))
                elif var in quantites:
                    varlist_dict[pipeline].append(var)
        if len(varlist_dict[pipeline]) != len(default_variables):
            print("Requried variables: {}".format(default_variables))
            print("Found variables: {}".format(varlist_dict))
            raise Exception("Not all variables found in tree, aborting")
        # Check availability of required input variables
        if category_mode != "inclusive":
            varlist_dict[pipeline].append("%s_max_index" % channel)
        if expression not in config['fraction_binning'].keys():
            varlist_dict[pipeline].append(expression)

        # one fake factor per tau is needed
        suffix_dict = {
            "et": [2],
            "mt": [2],
            "tt": [1, 2]
        }
        # for shifts, append the shift name to the suffix
        shift_suffix = ""
        if pipeline != "nominal":
            shift_suffix = "__%s" % pipeline
        for suffix in suffix_dict[channel]:

            output_buffer["nom_{}{}".format(suffix, shift_suffix)] = numpy.zeros(1, dtype=float)
            output_buffer["onlyqcd_{}{}".format(suffix, shift_suffix)] = numpy.zeros(1, dtype=float)
            if channel in ["mt","et"]:
                output_buffer["onlyw_{}{}".format(suffix, shift_suffix)] = numpy.zeros(1, dtype=float)
                output_buffer["onlytt_{}{}".format(suffix, shift_suffix)] = numpy.zeros(1, dtype=float)
            output_buffer["fracw_{}{}".format(suffix, shift_suffix)] = numpy.zeros(1, dtype=float)
            output_buffer["fracqcd_{}{}".format(suffix, shift_suffix)] = numpy.zeros(1, dtype=float)
            output_buffer["fractt_{}{}".format(suffix, shift_suffix)] = numpy.zeros(1, dtype=float)

            output_tree.Branch("ff{}_nom{}".format(suffix, shift_suffix), output_buffer["nom_{}{}".format(suffix, shift_suffix)], "ff{}_nom{}/D".format(suffix, shift_suffix))
            output_tree.Branch("ff{}_onlyqcd{}".format(suffix, shift_suffix), output_buffer["onlyqcd_{}{}".format(suffix, shift_suffix)], "ff{}_onlyqcd{}/D".format(suffix, shift_suffix))
            if channel in ["mt","et"]:
                output_tree.Branch("ff{}_onlyw{}".format(suffix, shift_suffix), output_buffer["onlyw_{}{}".format(suffix, shift_suffix)], "ff{}_onlyw{}/D".format(suffix, shift_suffix))
                output_tree.Branch("ff{}_onlytt{}".format(suffix, shift_suffix), output_buffer["onlytt_{}{}".format(suffix, shift_suffix)], "ff{}_onlytt{}/D".format(suffix, shift_suffix))

            output_tree.Branch("ff{}_fracw{}".format(suffix, shift_suffix), output_buffer["fracw_{}{}".format(suffix, shift_suffix)], "ff{}_fracw{}/D".format(suffix, shift_suffix))
            output_tree.Branch("ff{}_fracqcd{}".format(suffix, shift_suffix), output_buffer["fracqcd_{}{}".format(suffix, shift_suffix)], "ff{}_fracqcd{}/D".format(suffix, shift_suffix))
            output_tree.Branch("ff{}_fractt{}".format(suffix, shift_suffix), output_buffer["fractt_{}{}".format(suffix, shift_suffix)], "ff{}_fractt{}/D".format(suffix, shift_suffix))

            if pipeline == "nominal":
                for syst in unc_shifts[channel]:
                    for shift in ["up", "down"]:
                        output_buffer["%s_%s_%i" % (syst, shift, suffix)] = numpy.zeros(1, dtype=float)
                        output_tree.Branch(
                            "ff%i_%s_%s" % (suffix, syst, shift),
                            output_buffer["%s_%s_%i" % (syst, shift, suffix)],
                            "ff%i_%s_%s/D" % (suffix, syst, shift)
                        )
    # for _pipeline in required_piplines:
    #     print("Pipeline: {}".format(_pipeline))
    #     print("Variables: {}".format(varlist_dict[_pipeline]))
    # Fill tree
    print("eventrange", eventrange)
    for evt_i, event in enumerate(input_tree):
        if eventrange is not None:
            if evt_i < eventrange[0]:
                continue
            elif evt_i > eventrange[1] and eventrange[1] >= 0:  # latter condition allows to set negative upper limit in order to have it ignored
                print("Done with event range {}".format(eventrange))
                break
        eventcounter = evt_i - eventrange[0] if eventrange is not None else evt_i
        if(eventcounter % 100 == 0):
            print("Processing event {} / {} ({} %)".format(1 + eventcounter, 1 + eventrange[1] - eventrange[0], 100 * eventcounter / (1 + eventrange[1] - eventrange[0])))
        for pipeline in required_piplines:
            varlist = varlist_dict[pipeline]
            variableobjects = {}
            for i, var in enumerate(default_variables):
                variableobjects[var] = getattr(event, varlist[i])
            # print("Pipeline: {}".format(pipeline))
            # print('Using variables: %s' % varlist)
            # print('Using defaults: %s' % default_variables)
            # print("Variableobjects: {}".format(variableobjects))
            # for shifts, append the shift name to the suffix

            shift_suffix = ""
            if pipeline != "nominal":
                shift_suffix = "__%s" % pipeline

            for suffix in suffix_dict[channel]:
                inputs = []
                if category_mode == "inclusive":
                    cat_index = -1
                else:
                    if channel == "tt" and suffix == 2:
                        cat_index = int(getattr(event, "%s_max_index" % channel) + (0.5 * len(categories[channel])))
                    else:
                        cat_index = int(getattr(event, "%s_max_index" % channel) + 0.0)
                qcd_fraction = 0.0
                w_fraction = 0.0
                tt_fraction = 0.0
                if use_fractions_from_workspace:
                    fractions.var("cat").setVal(-1 if category_mode == "inclusive" else getattr(event, "%s_max_index" % channel))
                    fractions.var("m_vis").setVal(variableobjects["m_vis"])
                    fractions.var("njets").setVal(variableobjects["njets"])
                    fractions.var("aiso").setVal(suffix)
                    qcd_fraction = fractions.function("%s_frac_qcd" % channel[0]).getVal()
                    w_fraction = fractions.function("%s_frac_w" % channel[0]).getVal()
                    tt_fraction = fractions.function("%s_frac_tt" % channel[0]).getVal()
                else:
                    varvalue = 0.0
                    if expression == "njets_mvis":
                        varvalue = 300.0 * min(variableobjects["njets"], 2.0) + min(290.0, variableobjects["m_vis"])
                    elif expression == "njets2bins_mt_1_puppi":
                        varvalue = 180.0 * min(variableobjects["njets"], 1.0) + min(160.0, variableobjects["mt_1"])
                    elif expression == "njets3bins_mt_1_puppi":
                        varvalue = 180.0 * min(variableobjects["njets"], 2.0) + min(160.0, variableobjects["mt_1"])
                    elif expression == "njets2bins_m_vis":
                        varvalue = 250.0 * min(variableobjects["njets"], 1.0) + min(240.0, variableobjects["m_vis"])
                    elif expression == "njets3bins_m_vis":
                        varvalue = 250.0 * min(variableobjects["njets"], 2.0) + min(240.0, variableobjects["m_vis"])
                    elif expression in config['fraction_binning'].keys():
                        varvalue = eval(config['fraction_binning'][expression][channel]['expression'], {'min': min, 'max': max, 'm_vis': variableobjects["m_vis"], 'njets': variableobjects["njets"], 'mt_1_puppi': variableobjects["mt_1"]})
                    else:
                        varvalue = getattr(event, expression)
                    cat_fractions = fractions[channel][categories[channel][cat_index]]
                    bin_index = cat_fractions["data"].GetXaxis().FindBin(varvalue)
                    qcd_fraction = cat_fractions["QCD"].GetBinContent(bin_index)
                    w_fraction = cat_fractions["W"].GetBinContent(bin_index)

                    if channel in ["mt","et"]:
                        w_fraction_error = sqrt(power(cat_fractions["W"].GetBinError(bin_index)/cat_fractions["W"].GetBinContent(bin_index),2)+power(0.05,2))
                    tt_fraction = cat_fractions["TT"].GetBinContent(bin_index)

                fraction_sum = qcd_fraction + w_fraction + tt_fraction
                if fraction_sum == 0:
                    logger.critical("Division by zero -> forcing fraction_sum to 1")
                    fraction_sum = 1
                if not use_fractions_from_workspace and channel in ["mt","et"]:
                    w_fraction_w_up = (1.+w_fraction_error)*w_fraction
                    w_fraction_w_down = (1.-w_fraction_error)*w_fraction
                    w_fraction_w_up /= fraction_sum
                    w_fraction_w_down /= fraction_sum
                qcd_fraction /= fraction_sum
                w_fraction /= fraction_sum
                tt_fraction /= fraction_sum

                if not use_fractions_from_workspace and channel in ["mt","et"]:
                    qcd_fraction_w_up = 1 - tt_fraction - w_fraction_w_up
                    qcd_fraction_w_down = 1 - tt_fraction - w_fraction_w_down

                if channel == "tt":
                    inputs = [
                        variableobjects["pt_{}".format(suffix)],
                        variableobjects["pt_{}".format(3 - suffix)],
                        variableobjects["decaymode_{}".format(suffix)],
                        variableobjects["njets"],
                        variableobjects["m_vis"],
                        qcd_fraction,
                        w_fraction,
                        tt_fraction
                    ]
                else:
                    inputs = [
                        variableobjects["pt_2"], variableobjects["deltaR_ditaupair"], variableobjects["njets"], variableobjects["m_vis"], variableobjects["pt_1"],
                        variableobjects["mt_1"], variableobjects["iso_1"],
                        qcd_fraction,
                        w_fraction,
                        tt_fraction,
                        w_fraction_w_up, w_fraction_w_down, qcd_fraction_w_up, qcd_fraction_w_down
                    ]
                    if use_mva_tauid:
                        inputs = [
                            variableobjects["pt_2"], variableobjects["decaymode_2"],
                            variableobjects["njets"], variableobjects["m_vis"],
                            variableobjects["mt_1"], variableobjects["iso_1"],
                            qcd_fraction,
                            w_fraction,
                            tt_fraction
                        ]

                output_buffer["nom_{}{}".format(suffix, shift_suffix)][0] = ff.value(len(inputs), array('d', inputs))
                output_buffer["onlyqcd_{}{}".format(suffix, shift_suffix)][0] = ff.value(len(inputs), array('d', inputs), "ff_onlyqcd")
                if channel in ["mt","et"]:
                    output_buffer["onlyw_{}{}".format(suffix, shift_suffix)][0] = ff.value(len(inputs), array('d', inputs), "ff_onlyw")
                    output_buffer["onlytt_{}{}".format(suffix, shift_suffix)][0] = ff.value(len(inputs), array('d', inputs), "ff_onlytt")
                output_buffer["fracw_{}{}".format(suffix, shift_suffix)][0] = ff.value(len(inputs), array('d', inputs), "ff_fracw")
                output_buffer["fracqcd_{}{}".format(suffix, shift_suffix)][0] = ff.value(len(inputs), array('d', inputs), "ff_fracqcd")
                output_buffer["fractt_{}{}".format(suffix, shift_suffix)][0] = ff.value(len(inputs), array('d', inputs), "ff_fractt")

                if not (output_buffer["nom_{}{}".format(suffix, shift_suffix)][0]  >= 0.0 and output_buffer["nom_{}{}".format(suffix, shift_suffix)][0]  <= 999.0):
                    logger.info("Got invalid nominal weight %s for inputs [%s]" % (str(output_buffer["nom_{}{}".format(suffix, shift_suffix)][0] ), ', '.join(str(e) for e in inputs)))
                    output_buffer["nom_{}{}".format(suffix, shift_suffix)][0]  = 0.0
                if pipeline == "nominal":
                    for syst in unc_shifts[channel]:
                        for shift in ["up", "down"]:
                            output_buffer["%s_%s_%i" % (syst, shift, suffix)][0] = ff.value(len(inputs), array('d', inputs), "%s_%s" % (syst, shift))
                            if not (output_buffer["%s_%s_%i" % (syst, shift, suffix)][0] >= 0.0 and output_buffer["%s_%s_%i" % (syst, shift, suffix)][0] <= 999.0):
                                logger.info("Got invalid weight %s for syst shift %s and inputs [%s]" % (str(output_buffer["%s_%s_%i" % (syst, shift, suffix)][0]), syst, ', '.join(str(e) for e in inputs)))
                                output_buffer["%s_%s_%i" % (syst, shift, suffix)][0] = 0.0
            # ff.Delete() ; print 'bye'; exit(1)
        output_tree.Fill()

    # Save
    output_tree.Write()
    print("Successfully finished %s" % outputfile)

    # Clean up
    ff.Delete()
    ff_file.Close()
    # for input_friend_file in friendfiles:
    #     input_friend_file.Close()
    input_file.Close()
    output_file.Close()
    print('done')
    return 0



def load_fractions(configpath, configkey, use_fractions_from_workspace, workspace, fractionfiles, era):
    if use_fractions_from_workspace:
        logger.info("Loading workspace from %s" % workspace)
        f = ROOT.TFile(workspace)
        fractions = f.Get("w")
        f.Close()
        return fractions
    else:
        return determine_fractions(configpath, configkey, fractionfiles, era)

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
                logger.debug("Determining fractions: #{ch}#{ch}_{cat}#data_obs#{tag}#{era}#{expr}#125#".format(
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
                            logger.debug("Rescaled %s fraction to %f" % (fraction, subsubdict[fraction].GetBinContent(i)))
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
            "ff_qcd_syst", "ff_qcd_dm0_njet0_stat", "ff_qcd_dm0_njet1_stat", "ff_qcd_dm0_njet2_stat",
            "ff_w_syst", "ff_w_dm0_njet0_stat", "ff_w_dm0_njet1_stat", "ff_w_dm0_njet2_stat",
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
            "ff_qcd_dm0_njet0_morphed_stat", "ff_qcd_dm0_njet1_morphed_stat", "ff_qcd_dm0_njet2_morphed_stat",
            "ff_w_dm0_njet0_morphed_stat", "ff_w_dm0_njet1_morphed_stat", "ff_w_dm0_njet2_morphed_stat",
            "ff_tt_dm0_njet0_morphed_stat", "ff_tt_dm0_njet1_morphed_stat",
        ],
        "mt": [
            "ff_qcd_syst", "ff_qcd_dm0_njet0_stat", "ff_qcd_dm0_njet1_stat", "ff_qcd_dm0_njet2_stat",
            "ff_w_syst", "ff_w_dm0_njet0_stat", "ff_w_dm0_njet1_stat", "ff_w_dm0_njet2_stat",
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
            "ff_qcd_dm0_njet0_morphed_stat", "ff_qcd_dm0_njet1_morphed_stat", "ff_qcd_dm0_njet2_morphed_stat",
            "ff_w_dm0_njet0_morphed_stat", "ff_w_dm0_njet1_morphed_stat", "ff_w_dm0_njet2_morphed_stat",
            "ff_tt_dm0_njet0_morphed_stat", "ff_tt_dm0_njet1_morphed_stat", 
        ],
        "tt": [
            "ff_qcd_syst", "ff_qcd_dm0_njet0_stat", "ff_qcd_dm0_njet1_stat", "ff_qcd_dm0_njet2_stat",
            "ff_w_syst", "ff_tt_syst", "ff_w_frac_syst", "ff_tt_frac_syst",
            "ff_qcd_mvis",
            "ff_corr_qcd_mvis_osss",
            "ff_qcd_mvis_osss",
            "ff_qcd_tau2_pt",
            "ff_corr_qcd_mvis",
            "ff_corr_qcd_tau2_pt",
            "ff_qcd_mc",
            "ff_qcd_dm0_njet0_morphed_stat", "ff_qcd_dm0_njet1_morphed_stat", "ff_qcd_dm0_njet2_morphed_stat",
        ]
    }

    # Determine channel
    channels = ["et", "mt", "tt"]
    pipelines = ["nominal"]
    filename = os.path.basename(datafile)

    if "SingleElectron" in filename or "_ElTau" in filename:
        channels = ["et"]
    elif "SingleMuon" in filename or "_MuTau" in filename:
        channels = ["mt"]
    elif ("Tau" in filename and "Run%s" % era in filename) or "_TauTau" in filename:
        channels = ["tt"]
    else:
        channels = ["et", "mt", "tt"]
    if "Run%s" % era not in filename:
        pipelines = ["nominal", "tauEsOneProngUp", "tauEsOneProngDown",
                     "tauEsOneProngOnePiZeroUp", "tauEsOneProngOnePiZeroDown",
                     "tauEsThreeProngUp", "tauEsThreeProngDown",
                     "tauEsThreeProngOnePiZeroUp", "tauEsThreeProngOnePiZeroDown"]
    for channel in channels:
        for pipeline in pipelines:

            if pipeline_selection and pipeline_selection != "%s_%s" % (channel, pipeline):
                logger.debug('skip: %s because %s != %s_%s' % (pipeline, pipeline_selection, channel, pipeline))
                continue
            logger.debug('Processing pipeline: %s' % pipeline)

            # Prepare data inputs
            input_file = ROOT.TFile(datafile, "READ")
            input_tree = input_file.Get("%s_%s/%s" % (channel, pipeline, treename))
            friendfiles = []
            for input_friend_file in friendfilelists[channel]:
                friendfiles.append(ROOT.TFile(input_friend_file, "READ"))
                input_friend = friendfiles[-1].Get("%s_%s/%s" % (channel, pipeline, treename))
                input_tree.AddFriend(input_friend)

            # Check availability of required input variables
            varlist = ["pt_1", "pt_2", "decayMode_1", "decayMode_2", "m_vis", "njets", "iso_1", "mt_1", "mt_1_puppi","jpt_1"]
            if category_mode != "inclusive":
                varlist.append("%s_max_index" % channel)
            if expression not in config['fraction_binning'].keys():
                varlist.append(expression)

            for variable in varlist:
                if not input_tree.GetBranchStatus(variable):
                    logger.critical("Tree does not contain input variable '%s'" % variable)
                    raise Exception

            # Load fake factors histograms
            ff_file = ROOT.TFile.Open(fakefactordirectories[channel])
            print fakefactordirectories[channel]
            ff = ff_file.Get('ff_comb')

            # Prepare output
            logger.debug("...initialize %s" % outputfile)
            if not os.path.exists(os.path.dirname(outputfile)):
                os.mkdir(os.path.dirname(outputfile))
            output_file = ROOT.TFile(outputfile, rootfilemode)
            output_root_dir = output_file.mkdir("%s_%s" % (channel, pipeline))
            output_root_dir.cd()
            output_tree = ROOT.TTree(treename, treename)

            # one fake factor per tau is needed
            suffix = {
                "et": [2],
                "mt": [2],
                "tt": [1, 2]
            }
            output_buffer = {}
            for x in suffix[channel]:
                output_buffer["nom_%i" % x] = numpy.zeros(1, dtype=float)
                output_tree.Branch("ff%i_nom" % x, output_buffer["nom_%i" % x], "ff%i_nom/D" % x)
                output_buffer["onlyqcd_%i" % x] = numpy.zeros(1, dtype=float)
                if channel in ["mt","et"]:
                    output_buffer["onlyw_%i" % x] = numpy.zeros(1, dtype=float)
                    output_buffer["onlytt_%i" % x] = numpy.zeros(1, dtype=float)
                output_buffer["fracw_%i" % x] = numpy.zeros(1, dtype=float)
                output_buffer["fracqcd_%i" % x] = numpy.zeros(1, dtype=float)
                output_buffer["fractt_%i" % x] = numpy.zeros(1, dtype=float)

                output_tree.Branch("ff%i_nom" % x, output_buffer["nom_%i" % x], "ff%i_nom/D" % x)
                output_tree.Branch("ff%i_onlyqcd" % x, output_buffer["onlyqcd_%i" % x], "ff%i_onlyqcd/D" % x)
                if channel in ["mt","et"]:
                    output_tree.Branch("ff%i_onlyw" % x, output_buffer["onlyw_%i" % x], "ff%i_onlyw/D" % x)
                    output_tree.Branch("ff%i_onlytt" % x, output_buffer["onlytt_%i" % x], "ff%i_onlytt/D" % x)

                output_tree.Branch("ff%i_fracw" % x, output_buffer["fracw_%i" % x], "ff%i_fracw/D" % x)
                output_tree.Branch("ff%i_fracqcd" % x, output_buffer["fracqcd_%i" % x], "ff%i_fracqcd/D" % x)
                output_tree.Branch("ff%i_fractt" % x, output_buffer["fractt_%i" % x], "ff%i_fractt/D" % x)


                for syst in unc_shifts[channel]:
                    for shift in ["up", "down"]:
                        output_buffer["%s_%s_%i" % (syst, shift, x)] = numpy.zeros(1, dtype=float)
                        output_tree.Branch(
                            "ff%i_%s_%s" % (x, syst, shift),
                            output_buffer["%s_%s_%i" % (syst, shift, x)],
                            "ff%i_%s_%s/D" % (x, syst, shift)
                        )

            # Fill tree
            for evt_i, event in enumerate(input_tree):
                if eventrange is not None:
                    if evt_i < eventrange[0]:
                        continue
                    elif evt_i > eventrange[1] and eventrange[1] >= 0:  # latter condition allows to set negative upper limit in order to have it ignored
                        break
                for x in suffix[channel]:
                    inputs = []
                    if category_mode == "inclusive":
                        cat_index = -1
                    else:
                        if channel == "tt" and x == 2:
                            cat_index = int(getattr(event, "%s_max_index" % channel) + (0.5 * len(categories[channel])))
                        else:
                            cat_index = int(getattr(event, "%s_max_index" % channel) + 0.0)
                    qcd_fraction = 0.0
                    w_fraction = 0.0
                    tt_fraction = 0.0
                    if use_fractions_from_workspace:
                        fractions.var("cat").setVal(-1 if category_mode == "inclusive" else getattr(event, "%s_max_index" % channel))
                        fractions.var("m_vis").setVal(event.m_vis)
                        fractions.var("njets").setVal(event.njets)
                        fractions.var("aiso").setVal(x)
                        qcd_fraction = fractions.function("%s_frac_qcd" % channel[0]).getVal()
                        w_fraction = fractions.function("%s_frac_w" % channel[0]).getVal()
                        tt_fraction = fractions.function("%s_frac_tt" % channel[0]).getVal()
                    else:
                        varvalue = 0.0
                        if expression == "njets_mvis":
                            varvalue = 300.0 * min(event.njets, 2.0) + min(290.0, event.m_vis)
                        elif expression == "njets2bins_mt_1_puppi":
                            varvalue = 180.0 * min(event.njets, 1.0) + min(160.0, event.mt_1_puppi)
                        elif expression == "njets3bins_mt_1_puppi":
                            varvalue = 180.0 * min(event.njets, 2.0) + min(160.0, event.mt_1_puppi)
                        elif expression == "njets2bins_m_vis":
                            varvalue = 250.0 * min(event.njets, 1.0) + min(240.0, event.m_vis)
                        elif expression == "njets3bins_m_vis":
                            varvalue = 250.0 * min(event.njets, 2.0) + min(240.0, event.m_vis)
                        elif expression == "nbtag3bins_mt_1_puppi":
                            varvalue = 180.0 * min(event.nbtag, 2.0) + min(160.0, event.mt_1_puppi)
                        elif expression == "nbtag3bins_m_vis":
                            varvalue = 250.0 * min(event.nbtag, 2.0) + min(240.0, event.m_vis)
                        elif expression in config['fraction_binning'].keys():
                            varvalue = eval(config['fraction_binning'][expression][channel]['expression'], {'min': min, 'max': max, 'm_vis': event.m_vis, 'njets': event.njets, 'mt_1_puppi': event.mt_1_puppi})
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
                            getattr(event, "pt_%i" % x),
                            getattr(event, "pt_%i" % (3 - x)),
                            getattr(event, "decayMode_%i" % x), event.njets,
                            event.m_vis, event.jpt_1, qcd_fraction,
                            w_fraction,
                            tt_fraction
                        ]
                    else:
                        inputs = [
                            event.pt_2, event.decayMode_2, event.njets, event.m_vis, event.pt_1,
                            event.mt_1, event.iso_1,
                            qcd_fraction,
                            w_fraction,
                            tt_fraction,
                            w_fraction_w_up, w_fraction_w_down, qcd_fraction_w_up, qcd_fraction_w_down
                        ]
                        if use_mva_tauid:
                            inputs = [
                                event.pt_2, event.decayMode_2,
                                event.njets, event.m_vis,
                                event.mt_1, event.iso_1,
                                qcd_fraction,
                                w_fraction,
                                tt_fraction
                            ]

                    output_buffer["nom_%i" % x][0] = ff.value(len(inputs), array('d', inputs))
                    output_buffer["onlyqcd_%i" % x][0] = ff.value(len(inputs), array('d', inputs), "ff_onlyqcd")
                    if channel in ["mt","et"]:
                        output_buffer["onlyw_%i" % x][0] = ff.value(len(inputs), array('d', inputs), "ff_onlyw")
                        output_buffer["onlytt_%i" % x][0] = ff.value(len(inputs), array('d', inputs), "ff_onlytt")
                    output_buffer["fracw_%i" % x][0] = ff.value(len(inputs), array('d', inputs), "ff_fracw")
                    output_buffer["fracqcd_%i" % x][0] = ff.value(len(inputs), array('d', inputs), "ff_fracqcd")
                    output_buffer["fractt_%i" % x][0] = ff.value(len(inputs), array('d', inputs), "ff_fractt")

                    if not (output_buffer["nom_%i" % x][0] >= 0.0 and output_buffer["nom_%i" % x][0] <= 999.0):
                        logger.info("Got invalid nominal weight %s for inputs [%s]" % (str(output_buffer["nom_%i" % x][0]), ', '.join(str(e) for e in inputs)))
                        output_buffer["nom_%i" % x][0] = 0.0
                    for syst in unc_shifts[channel]:
                        for shift in ["up", "down"]:
                            output_buffer["%s_%s_%i" % (syst, shift, x)][0] = ff.value(len(inputs), array('d', inputs), "%s_%s" % (syst, shift))
                            if not (output_buffer["%s_%s_%i" % (syst, shift, x)][0] >= 0.0 and output_buffer["%s_%s_%i" % (syst, shift, x)][0] <= 999.0):
                                logger.info("Got invalid weight %s for syst shift %s and inputs [%s]" % (str(output_buffer["%s_%s_%i" % (syst, shift, x)][0]), syst, ', '.join(str(e) for e in inputs)))
                                output_buffer["%s_%s_%i" % (syst, shift, x)][0] = 0.0
                # ff.Delete() ; print 'bye'; exit(1)
                output_tree.Fill()

            # Save
            output_tree.Write()
            logger.debug("Successfully finished %s" % outputfile)

            # Clean up
            ff.Delete()
            ff_file.Close()
            for input_friend_file in friendfiles:
                input_friend_file.Close()
            input_file.Close()
            output_file.Close()
            print 'done'
            return 0

    print 'done (no matching pipelines found)'
    return 1


def load_fractions(configpath, configkey, use_fractions_from_workspace, workspace, fractionfiles, era):
    if use_fractions_from_workspace:
        logger.info("Loading workspace from %s" % workspace)
        f = ROOT.TFile(workspace)
        fractions = f.Get("w")
        f.Close()
        return fractions
    else:
        return determine_fractions(configpath, configkey, fractionfiles, era)

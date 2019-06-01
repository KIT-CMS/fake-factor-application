#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import yaml
import ROOT
import numpy
import copy
from array import array

import logging
logger = logging.getLogger()


def determine_fractions(configpath, configkey, fractionsfile, era):
    config = yaml.load(open(configpath))
    if not "categories" in config.keys():
        logger.critical("Config file %s has to contain key 'categories'!" % configpath)
        raise Exception
    categories = config["categories"]
    hist_file = ROOT.TFile(fractionsfile, "READ")
    era_labels = {
        "2016" : "Run2016",
        "2017" : "Run2017ReReco31Mar"
            }
    composition = {
        "mt": {
            "data": ["data_obs"],
            "W": ["W", "VVJ", "ZJ"],
            "TT": ["TTJ"],
            "QCD": ["QCD"],
            "real": ["EMB", "ZL", "TTL", "VVL"] #["ZTT", "ZL", "TTT", "TTL", "VVT", "VVL"]
        },
        "et": {
            "data": ["data_obs"],
            "W": ["W", "VVJ", "ZJ"],
            "TT": ["TTJ"],
            "QCD": ["QCD"],
            "real": ["EMB", "ZL", "TTL", "VVL"] #["ZTT", "ZL", "TTT", "TTL", "VVT", "VVL"]
        },
        "tt": {
            "data": ["data_obs"],
            "W": ["W", "VVJ", "ZJ"],
            "TT": ["TTJ"],
            "QCD": ["QCD"],
            "real": ["EMB", "ZL", "TTL", "VVL"] #["ZTT", "ZL", "TTT", "TTL", "VVT", "VVL"]
        }
    }
    fractions = {}
    for channel in categories.keys():
        subdict = {}
        for category in categories[channel]:
            subsubdict = {}
            for fraction in composition[channel].keys():
                subsubdict[fraction] = copy.deepcopy(
                    hist_file.Get(
                        "#{ch}#{ch}_{cat}#data_obs#smhtt#{era}#{expr}#125#".
                        format(
                            ch=channel,
                            cat=category,
                            era=era_labels[era],
                            expr=configkey)))
                subsubdict[fraction].Reset()
            for fraction in composition[channel].keys():
                if fraction == "QCD":
                    continue
                for process in composition[channel][fraction]:
                    subsubdict[fraction].Add(
                        hist_file.Get(
                            "#{ch}#{ch}_{cat}#{proc}#smhtt#{era}#{expr}#125#".
                            format(
                                ch=channel,
                                cat=category,
                                proc=process,
                                era=era_labels[era],
                                expr=configkey)))
                if fraction == "data":
                    subsubdict["QCD"].Add(subsubdict[fraction], 1.0)
                else:
                    subsubdict["QCD"].Add(subsubdict[fraction], -1.0)
            #normalize to data to get fractions
            denominator_hist = copy.deepcopy(subsubdict["data"])
            for fraction in composition[channel].keys():
                subsubdict[fraction].Divide(denominator_hist)
            #sanity check: if QCD negative i.e. data < MC, normalize to MC
            for i in range(subsubdict["QCD"].GetNbinsX()+2):
                qcd_fraction = subsubdict["QCD"].GetBinContent(i)
                if qcd_fraction < 0.0:
                    logger.info("Found bin with negative QCD fraction (%s, %s, index %i). Set fraction to zero and rescale other fractions..."%(channel, category, i))
                    subsubdict["QCD"].SetBinContent(i, 0)
                    for fraction in composition[channel].keys():
                        if not fraction == "data":
                            subsubdict[fraction].SetBinContent(i, subsubdict[fraction].GetBinContent(i) / (1.0 - qcd_fraction))
                            logger.debug("Rescaled %s fraction to %f"%(fraction, subsubdict[fraction].GetBinContent(i)))
            subdict[category] = subsubdict
        fractions[channel] = subdict
    hist_file.Close()
    return fractions


def apply_fake_factors(datafile, friendfilelists, outputfile, category_mode, fakefactordirectories, fractions, use_fractions_from_worspace, configpath, expression, era, pipeline_selection=None, treename="ntuple", eventrange=None, rootfilemode="update"):
    
    config = yaml.load(open(configpath))
    if not "categories" in config.keys():
        logger.critical("Config file %s has to contain key 'categories'!" % configpath)
        raise Exception
    categories = config["categories"]

    rootfilemode = rootfilemode.lower()
    if rootfilemode not in ["update", "recreate"]:
        logger.critical("Mode %s not appropriate for create of ROOT file. Please choose from 'update' and 'recreate'"%rootfilemode)
        raise Exception

    unc_shifts = {  #documented in https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauJet2TauFakes
        "et": [
            "ff_qcd_syst", "ff_qcd_dm0_njet0_stat", "ff_qcd_dm0_njet1_stat",
            "ff_w_syst", "ff_w_dm0_njet0_stat", "ff_w_dm0_njet1_stat",
            "ff_tt_syst", "ff_tt_dm0_njet0_stat", "ff_tt_dm0_njet1_stat"
        ],
        "mt": [
            "ff_qcd_syst", "ff_qcd_dm0_njet0_stat", "ff_qcd_dm0_njet1_stat",
            "ff_w_syst", "ff_w_dm0_njet0_stat", "ff_w_dm0_njet1_stat",
            "ff_tt_syst", "ff_tt_dm0_njet0_stat", "ff_tt_dm0_njet1_stat"
        ],
        "tt": [
            "ff_qcd_syst", "ff_qcd_dm0_njet0_stat", "ff_qcd_dm0_njet1_stat",
            "ff_w_syst", "ff_tt_syst", "ff_w_frac_syst", "ff_tt_frac_syst"
        ]
    }

    #determine channel
    channels = ["et", "mt", "tt"]
    pipelines = ["nominal"]
    filename=os.path.basename(datafile)
    if "SingleElectron" in filename or "_ElTau" in filename:
        channels = ["et"]
    elif "SingleMuon" in filename or "_MuTau" in filename:
        channels = ["mt"]
    elif ("Tau" in filename and "Run%s"%era in filename) or "_TauTau" in filename:
        channels = ["tt"]
    else:
        channels = ["et", "mt", "tt"]
    if not "Run%s"%era in filename:
        pizero="PiZeros" if era=="2016" else "OnePiZero"
        pipelines = ["nominal", "tauEsOneProngUp", "tauEsOneProngDown",
                     "tauEsOneProng%sUp"%pizero, "tauEsOneProng%sDown"%pizero,
                     "tauEsThreeProngUp", "tauEsThreeProngDown"]
    for channel in channels:
        for pipeline in pipelines:
            if pipeline_selection and pipeline_selection!="%s_%s"%(channel, pipeline):
                continue
            
            #prepare data inputs
            input_file = ROOT.TFile(datafile, "READ")
            input_tree = input_file.Get("%s_%s/%s" %(channel, pipeline, treename))
            friendfiles = []
            for input_friend_file in friendfilelists[channel]:
                friendfiles.append(ROOT.TFile(input_friend_file, "READ"))
                input_friend = friendfiles[-1].Get("%s_%s/%s" %(channel, pipeline, treename))
                input_tree.AddFriend(input_friend)
            
            #check availability of required input variables
            varlist =  ["pt_1", "pt_2", "decayMode_1", "decayMode_2", "m_vis", "njets", "iso_1", "mt_1"]
            if category_mode != "inclusive":
                varlist.append("%s_max_index" % channel)
            if expression != "njets_mvis":
                varlist.append(expression)
            for variable in varlist:
                if not input_tree.GetBranchStatus(variable):
                    logger.critical("Tree does not contain input variable '%s'" % variable)
                    raise Exception

            #load fake factors histograms
            ff_file = ROOT.TFile.Open(fakefactordirectories[channel])
            ff = ff_file.Get('ff_comb')

            #prepare output
            logger.debug(
                "...initialize %s" % outputfile)
            if not os.path.exists(os.path.dirname(outputfile)):
                os.mkdir(os.path.dirname(outputfile))
            output_file = ROOT.TFile(outputfile, rootfilemode)
            output_root_dir = output_file.mkdir("%s_%s" %(channel, pipeline))
            output_root_dir.cd()
            output_tree = ROOT.TTree(treename, treename)

            suffix = {  #one fake factor per tau is needed
                "et": [2],
                "mt": [2],
                "tt": [1, 2]
            }
            output_buffer = {}
            for x in suffix[channel]:
                output_buffer["nom_%i" % x] = numpy.zeros(1, dtype=float)
                output_tree.Branch("ff%i_nom" % x, output_buffer["nom_%i" % x],
                                "ff%i_nom/D" % x)
                for syst in unc_shifts[channel]:
                    for shift in ["up", "down"]:
                        output_buffer["%s_%s_%i" % (syst, shift, x)] = numpy.zeros(
                            1, dtype=float)
                        output_tree.Branch("ff%i_%s_%s" % (x, syst, shift),
                                        output_buffer["%s_%s_%i" %
                                                        (syst, shift,
                                                        x)], "ff%i_%s_%s/D" % (x, syst, shift))

            #fill tree
            for evt_i, event in enumerate(input_tree):
                if eventrange!=None:
                    if evt_i<eventrange[0]:
                        continue
                    elif evt_i>eventrange[1] and eventrange[1]>=0: #latter condition allows to set negative upper limit in order to have it ignored
                        break
                for x in suffix[channel]:
                    inputs = []
                    cat_index = -1 if category_mode=="inclusive" else int(
                        getattr(event, "%s_max_index" % channel) +
                        (0.5 * len(categories[channel])
                        if channel == "tt" and x == 2 else 0.0))
                    qcd_fraction = 0.0
                    w_fraction = 0.0            
                    tt_fraction = 0.0
                    if use_fractions_from_worspace:
                        fractions.var("cat").setVal(-1 if category_mode=="inclusive" else getattr(event, "%s_max_index" % channel))
                        fractions.var("m_vis").setVal(event.m_vis)
                        fractions.var("njets").setVal(event.njets)
                        fractions.var("aiso").setVal(x)
                        qcd_fraction = fractions.function("%s_frac_qcd"%channel[0]).getVal()
                        w_fraction = fractions.function("%s_frac_w"%channel[0]).getVal()
                        tt_fraction = fractions.function("%s_frac_tt"%channel[0]).getVal()
                    else:
                        varvalue = 0.0
                        if expression=="njets_mvis":
                            varvalue = 300.0*min(event.njets, 2.0) + min(290.0, event.m_vis)
                        else:
                            varvalue = getattr(event, expression)
                        cat_fractions = fractions[channel][categories[channel][cat_index]]
                        bin_index = cat_fractions["data"].GetXaxis().FindBin(varvalue)
                        qcd_fraction = cat_fractions["QCD"].GetBinContent(bin_index)
                        w_fraction = cat_fractions["W"].GetBinContent(bin_index)
                        tt_fraction = cat_fractions["TT"].GetBinContent(bin_index)
                    fraction_sum = qcd_fraction + w_fraction + tt_fraction
                    qcd_fraction /= fraction_sum
                    w_fraction /= fraction_sum
                    tt_fraction /= fraction_sum
                    if channel == "tt":
                        inputs = [
                            getattr(event, "pt_%i" % x),
                            getattr(event, "pt_%i" % (3 - x)),
                            getattr(event, "decayMode_%i" % x), event.njets,
                            event.m_vis, qcd_fraction,
                            w_fraction,
                            tt_fraction
                        ]
                    else:
                        inputs = [
                            event.pt_2, event.decayMode_2, event.njets, event.m_vis,
                            event.mt_1, event.iso_1,
                            qcd_fraction,
                            w_fraction,
                            tt_fraction
                        ]
                    output_buffer["nom_%i" % x][0] = ff.value(
                        len(inputs), array('d', inputs))
                    if not (output_buffer["nom_%i" % x][0] >= 0.0 and output_buffer["nom_%i" % x][0] <= 999.0):
                        logger.info("Got invalid nominal weight %s for inputs [%s]"%(str(output_buffer["nom_%i" % x][0]), ', '.join(str(e) for e in inputs)))
                        output_buffer["nom_%i" % x][0] = 0.0
                    for syst in unc_shifts[channel]:
                        for shift in ["up", "down"]:
                            output_buffer["%s_%s_%i" % (syst, shift, x)][0] = ff.value(
                                len(inputs), array('d', inputs), "%s_%s" % (syst,
                                                                            shift))
                            if not (output_buffer["%s_%s_%i" % (syst, shift, x)][0] >= 0.0 and output_buffer["%s_%s_%i" % (syst, shift, x)][0] <= 999.0):
                                logger.info("Got invalid weight %s for syst shift %s and inputs [%s]"%(str(output_buffer["%s_%s_%i" % (syst, shift, x)][0]), syst, ', '.join(str(e) for e in inputs)))
                                output_buffer["%s_%s_%i" % (syst, shift, x)][0] = 0.0
                output_tree.Fill()

            #save
            output_tree.Write()
            logger.debug("Successfully finished %s" % outputfile)

            #clean up
            ff.Delete()
            ff_file.Close()
            for input_friend_file in friendfiles:
                input_friend_file.Close()
            input_file.Close()
            output_file.Close()


def load_fractions(configpath, configkey, use_fractions_from_worspace, workspace, fractionfiles, era):
    if use_fractions_from_worspace:
        logger.info("Loading workspace from %s"%workspace)
        f = ROOT.TFile(workspace)
        fractions = f.Get("w")
        f.Close()
        return fractions
    else:
        return determine_fractions(configpath, configkey, fractionfiles, era)

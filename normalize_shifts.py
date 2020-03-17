#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import argparse
import logging
import copy
logger = logging.getLogger()


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Calculate fake factors and create friend trees.")
    parser.add_argument("input", type=str, help="Shape-producer output file.")
    return parser.parse_args()


def setup_logging(output_file, level=logging.DEBUG):
    logger.setLevel(level)
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    file_handler = logging.FileHandler(output_file, "w")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)


def main(args):
    file_ = ROOT.TFile(args.input, "UPDATE")
    for key in [k for k in file_.GetListOfKeys() if not "output_tree" == k.GetName()]:
        # Parse name
        name = key.GetName()

        if name[-1] == "#":  # shape is nominal shape
            continue

        split = [x for x in name.split("#") if not x == ""]
        channel = split[0]
        category = split[1]
        process = split[2]

        if not process == "jetFakes":  # ff uncertainties apply only on jetFakes process
            continue

        if "frac_w" in name or "_mc" in name or "tt_sf" in name: # true systematic uncertainties are not altered
            continue

        shift = split[7]
        if shift[-4:] == "Down":
            shift_type = "down"
            continue #run down together with up
        elif shift[-2:] == "Up":
            shift_type = "up"
        else:
            logger.critical("Cannot determine shift type of systematic %s.",
                            name)
            raise Exception

        # Renormalize if systematic has sub-string _ff_	
        if "_ff_" in name:
            h_shift = file_.Get(name)
            h_shift_down = file_.Get(name.replace("Up", "Down"))
            if h_shift == None or h_shift_down == None:
                logger.critical("Failed to get shape syst. histogram %s.",
                                name)
                raise Exception

            nominal = "#" + "#".join(split[:-1]) + "#"
            h_nominal = file_.Get(nominal)
            if h_nominal == None:
                logger.critical("Failed to get nominal histogram %s.", nominal)
                raise Exception

            norm_shift = h_shift.Integral()
            norm_shift_down = h_shift_down.Integral()
            norm_nominal = h_nominal.Integral()
            if norm_shift == 0 or norm_shift_down == 0:
                logger.warning("Found shift with integral of zero for systematic %s. Continue.",
                        name)
                continue
            scale = norm_nominal / norm_shift
            scale_down = norm_nominal / norm_shift_down
            logger.debug(
                "Renormalize systematic %s (%f) with integral of %s (%f): %f",
                name, norm_shift, nominal, norm_nominal, scale)
            h_shift.Scale(scale)
            h_shift.Write()
            logger.debug(
                "Renormalize systematic %s (%f) with integral of %s (%f): %f",
                name.replace("Up", "Down"), norm_shift, nominal, norm_nominal, scale)
            h_shift_down.Scale(scale_down)
            h_shift_down.Write()
    logger.info("Successfully normalized fake factor systematics to nominal.")
    file_.Close()


if __name__ == "__main__":
    args = parse_arguments()
    setup_logging("normalize_shifts.log", logging.INFO)
    main(args)

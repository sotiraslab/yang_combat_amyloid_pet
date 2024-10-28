#!/usr/bin/env python
# ==============================================================================

# Name: 3_harmonize.py
# Author: Braden Yang
# Created: 12/21/2023
# Description: Harmonize simulated SUVR data using Centiloid and ComBat

# ==============================================================================

# *** toggle interactive mode ***
INTERACTIVE = False
# *** toggle interactive mode ***

# ==========================
# ===== IMPORT MODULES =====
# ==========================

import argparse
from os import path, makedirs, chdir
import sys

import joblib
import numpy as np
import pandas as pd

# ======================
# ===== PARSE ARGS =====
# ======================

parser = argparse.ArgumentParser(
    description="Harmonize simulated SUVR data using Centiloid and ComBat",
    formatter_class=argparse.RawDescriptionHelpFormatter
)

parser.add_argument("-w", "--wdir", help="path to root directory of project repo")
parser.add_argument("-i", "--input", help="filepath of input CSV file with simulated data")
parser.add_argument("-o", "--odir", help="output directory")
parser.add_argument("-c", "--combat_dir", help="path to directory containing trained ComBat models as .joblib")
parser.add_argument("-f", "--feature_list", help="path to list of features to harmonize in ComBat model")

if INTERACTIVE:
    args = parser.parse_args(args = [])
    # define preset arguments
    args.input = "simulation/sim.csv"
    args.odir = "simulation"
else:
    args = parser.parse_args()

chdir(args.wdir)

if not path.exists(args.odir): makedirs(args.odir, exist_ok = True)  # make output directory

# =========================
# ===== IMPORT COMBAT =====
# =========================

sys.path.insert(0, path.join(args.wdir, "code/1_Harmonize"))
from combat import *

# ============================
# ===== DEFINE FUNCTIONS =====
# ============================

def convert_to_centiloid(suvr, tracer, study, pvc = False):

    # TODO: move to byyfunctions

    if study == "ADNI":

        if tracer == "AV45":
            cl = 196.9 * suvr - 196.03
        elif tracer == "FBB":
            cl = 159.08 * suvr - 151.65
        else:
            cl = None
    
    elif study == "OASIS":

        if pvc:
            if tracer == "AV45":
                cl = 53.6 * suvr - 43.2
            elif tracer == "PIB":
                cl = 45.0 * suvr - 47.5
            else:
                cl = None
        
        else:
            if tracer == "AV45":
                cl = 163.6 * suvr - 181.0
            elif tracer == "PIB":
                cl = 111.8 * suvr - 119.3
            else:
                cl = None
    
    else:
        cl = None

    return cl

def map_centiloid(X, feature_cols, tracer_col, study_col, pvc_col = None):

    # TODO: move to byyfunctions

    X_centiloid = X.copy()

    if pvc_col is None:

        # if no PVC col provided, assume non-PVC
        pvc_col = np.zeros((len(tracer_col),)).astype(bool)

    for f in feature_cols:

        X_centiloid.loc[:, f] = list(map(
            convert_to_centiloid,
            X.loc[:, f],
            tracer_col,
            study_col,
            pvc_col
        ))

    return X_centiloid  

# ============================
# ===== DEFINE VARIABLES =====
# ============================

# load list of features to harmonize
features_to_harmonize = np.loadtxt(args.feature_list, dtype = str)

# ===============================
# ===== LOAD SIMULATED DATA =====
# ===============================

sim_df = pd.read_csv(args.input)

# make dummy variables
sim_df = PETDataLoader._create_dummy_variables(
    sim_df,
    ["sex", "apoe", "trial.group.fct"]
)

# ================================
# ===== HARMONIZE: CENTILOID =====
# ================================

centiloid_df = map_centiloid(
    sim_df,
    features_to_harmonize,
    sim_df["tracer"],
    np.repeat("OASIS", sim_df.shape[0])
)
centiloid_df["harmonization_method"] = "cl"

# =============================
# ===== HARMONIZE: COMBAT =====
# =============================

# no covariates
combat_model__nocovar = joblib.load(path.join(args.combat_dir, "combat__no_covar.joblib"))
combat_nocovar_df = combat_model__nocovar.transform(sim_df)
combat_nocovar_df["harmonization_method"] = "combat_nocovar"

# age, sex, apoe
combat_model__age_sex_apoe = joblib.load(path.join(args.combat_dir, "combat__age_sex_apoe.joblib"))
combat_age_sex_apoe_df = combat_model__age_sex_apoe.transform(sim_df)
combat_age_sex_apoe_df["harmonization_method"] = "combat_age_sex_apoe"

# merge
sim_df["harmonization_method"] = "suvr"
harmonized_df = pd.concat(
    [sim_df, centiloid_df, combat_nocovar_df, combat_age_sex_apoe_df],
    axis = 0
)

# =====================
# ===== SAVE DATA =====
# =====================

harmonized_df.to_csv(
    path.join(args.odir, "sim_harm_cl_combat.csv"),
    sep = ",",
    index = False
)

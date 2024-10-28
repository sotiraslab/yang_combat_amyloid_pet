# ==============================================================================

# Name: 2_train_combat.py
# Author: Braden Yang
# Created: 12/21/2023
# Description: Train ComBat on amyloid-positive subject data for simulation
#   experiment

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
    description="Train ComBat on amyloid-positive subject data for simulation experiment",
    formatter_class=argparse.RawDescriptionHelpFormatter
)

parser.add_argument("-w", "--wdir", help="path to root directory of project repo")
parser.add_argument("-i", "--input", help="filepath of CSV containing amyloid-positive data")
parser.add_argument("-o", "--odir", help="path to output directory")
parser.add_argument("-f", "--feature_list", help="path to list of features to harmonize in ComBat model")

if INTERACTIVE:
    args = parser.parse_args(args = [])
    # define preset arguments
    args.input = "data/simulation/amypos.csv"
    args.odir = "data/simulation"
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
# ===== DEFINE VARIABLES =====
# ============================

# covariates
combat_covar = [
    "age",
    "sex_M",
    "apoe_1"
]

# load list of features to harmonize
features_to_harmonize = np.loadtxt(args.feature_list, dtype = str)

# ===========================
# ===== PRETRAIN COMBAT =====
# ===========================

train_df = pd.read_csv(args.input)
train_df["apoe"] = train_df["apoe"].astype("Int64")
train_df = PETDataLoader._create_dummy_variables(
    train_df,
    ["sex", "apoe"]
)

# train ComBat (no covariates)
combat_model_nocovar = GAMComBat(
    features = features_to_harmonize,
    covariates = None,
    smooth_terms = None,
    site_name = "tracer",
    random_state = 42
)
combat_model_nocovar.fit(train_df)
joblib.dump(combat_model_nocovar, path.join(args.odir, "combat__no_covar.joblib"))

# train ComBat with covariates
combat_model_covar = GAMComBat(
    features = features_to_harmonize,
    covariates = combat_covar,
    smooth_terms = None,
    site_name = "tracer",
    random_state = 42
)
combat_model_covar.fit(train_df)
joblib.dump(combat_model_covar, path.join(args.odir, f"combat__age_sex_apoe.joblib"))

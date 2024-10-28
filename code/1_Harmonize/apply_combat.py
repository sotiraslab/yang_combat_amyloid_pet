#!/usr/bin/env python
# ==============================================================================

# Name: apply_combat.py
# Author: Braden Yang
# Created: 10/04/2023
# Description: Train and apply ComBat on OASIS head-to-head crossover dataset

# ==============================================================================

# *** toggle interactive mode ***
INTERACTIVE = False
# *** toggle interactive mode ***

# ==========================
# ===== IMPORT MODULES =====
# ==========================

import argparse
from os import path, makedirs

import joblib

from combat import *

# ======================
# ===== PARSE ARGS =====
# ======================

parser = argparse.ArgumentParser(
    description="Train and apply ComBat on OASIS head-to-head crossover dataset",
    formatter_class=argparse.RawDescriptionHelpFormatter
)

parser.add_argument("-d", "--data_dir", help="path to data directory")
parser.add_argument("-o", "--odir", help="path to output directory")
parser.add_argument("-b", "--batch_other", default=None, help="name of additional batch effect to correct for; if not specified, then only tracer will be used as the batch effect")
parser.add_argument("-m", "--white_matter", action="store_true", default=False, help="add WM regions as features to ComBat model")

if INTERACTIVE:
    args = parser.parse_args(args = [])
    # define preset arguments
    args.batch_other = None
    args.white_matter = True
else:
    args = parser.parse_args()

if not path.exists(args.odir): makedirs(args.odir, exist_ok = True)  # make output directory

# =================================
# ===== LOAD AND PREPARE DATA =====
# =================================

# load data
pet_data = PETDataLoader(path.join(args.data_dir))
# pet_data.load_data()
pet_data.load_data(all_scans = True, batch_other = args.batch_other, white_matter = args.white_matter)
crossover_df = pet_data.crossover
crossover_train_df = pet_data.crossover_train

# feature names
feature_names = pet_data.feature_names
if args.white_matter:
    feature_names = np.concatenate([feature_names, pet_data.wm_feature_names])

# define covariates
covariates = ["age", "sex_M", "apoe_1"]

# filter nan SITE
crossover_df = crossover_df[~crossover_df["SITE"].isna()]
crossover_train_df = crossover_train_df[~crossover_train_df["SITE"].isna()]

# for pet_scanner_model only: filter out "1094" scanner
if args.batch_other == "pet_scanner_model":
    crossover_df = crossover_df[crossover_df["pet_scanner_model"] != "1094"]
    crossover_train_df = crossover_train_df[crossover_train_df["pet_scanner_model"] != "1094"]

# save data df for generating demographics table
demographics_dir = path.join(args.data_dir, "demographics")
if not path.exists(demographics_dir): makedirs(demographics_dir, exist_ok = True)

crossover_csv_path = path.join(demographics_dir, "crossover.csv")
if not path.exists(crossover_csv_path):
    crossover_df.to_csv(crossover_csv_path, index = False)

crossover_train_csv_path = path.join(demographics_dir, "crossover_train.csv")
if not path.exists(crossover_train_csv_path):
    crossover_train_df.to_csv(crossover_train_csv_path, index = False)

# ==================
# ===== COMBAT =====
# ==================

# ComBat with no covariates
combat_nocovar_model = GAMComBat(
    features = feature_names,
    covariates = None,
    smooth_terms = None,
    random_state = 42
)
combat_nocovar_model.fit(crossover_train_df)
combat_nocovar_df = combat_nocovar_model.transform(crossover_df)

# ComBat with linear age, sex, APOE
combat_linear_model = GAMComBat(
    features = feature_names,
    covariates = covariates,
    smooth_terms = None,
    smooth_bounds = None,
    random_state = 42
)
combat_linear_model.fit(crossover_train_df)
combat_linear_df = combat_linear_model.transform(crossover_df)

# GAM-ComBat with non-linear age, sex, APOE
combat_gam_model = GAMComBat(
    features = feature_names,
    covariates = covariates,
    smooth_terms = ["age"],
    smooth_bounds = get_smooth_bounds([crossover_df, crossover_train_df], "age"),
    random_state = 42
)
combat_gam_model.fit(crossover_train_df)
combat_gam_df = combat_gam_model.transform(crossover_df)

# REVISIONS: remove WM regions from results
if args.white_matter:
    combat_nocovar_df.drop(columns=pet_data.wm_feature_names, inplace=True)
    combat_linear_df.drop(columns=pet_data.wm_feature_names, inplace=True)
    combat_gam_df.drop(columns=pet_data.wm_feature_names, inplace=True)

# ========================
# ===== SAVE OUTPUTS =====
# ========================

def save_combat(df, model, prefix):
    
    df.to_csv(prefix + ".csv", index = False)
    dict2json(model.get_param_dict(), prefix + "_params.json")
    joblib.dump(model, prefix + "_model.joblib")

save_combat(combat_nocovar_df, combat_nocovar_model, path.join(args.odir, "combat_nocovar"))
save_combat(combat_linear_df, combat_linear_model, path.join(args.odir, "combat_linear"))
save_combat(combat_gam_df, combat_gam_model, path.join(args.odir, "combat_gam"))

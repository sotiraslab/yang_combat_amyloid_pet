
import json
from os import path

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedGroupKFold

import neuroHarmonize

# ======== classes ========

class GAMComBat(BaseEstimator, TransformerMixin):

    # https://www.andrewvillazon.com/custom-scikit-learn-transformers/

    def __init__(
        self,
        features,
        covariates,
        site_name = "SITE",
        smooth_terms = [],
        smooth_bounds = None,
        random_state = None,
        verbose = False
    ):
        """
        Parameters
        ----------
        features : list
            list of feature names in data df
        covariates : list
            list of covariate names in data df
        site_name : str
            name of column that encodes "SITE" information
        smooth_terms : list of str
            list of the smooth term(s) to use (if only one, still enlist)
        smooth_bounds : list of size (2,)
            bounds for smooth term
        random_state : int
            optional random state
        verbose : bool
            if True, print status statements
        """

        self.features = features
        self.covariates = covariates
        self.site_name = site_name
        self.smooth_terms = smooth_terms
        self.smooth_bounds = smooth_bounds
        self.random_state = random_state
        self.verbose = verbose

    def _vprint(self, s):

        if self.verbose: print(s)

    def _prepare_data(self, X, covar_only = False):

        # create new site column
        if not self.site_name == "SITE":
            X.loc[:, "SITE"] = X.loc[:, self.site_name]

        # separate into data and covariates
        if self.covariates is None:
            covar = X.loc[:, ["SITE"]]
        else:
            if not type(self.covariates) in [np.ndarray, list]:
                cov_list = [self.covariates]
            else:
                cov_list = self.covariates
            covar = X.loc[:, np.concatenate([cov_list, ["SITE"]])]

        if covar_only: return covar
        
        data = X.loc[:, self.features]

        return data, covar

    def fit(self, X, y = None):

        self._vprint("++ Training ComBat-GAM model on training data ++")

        # separate into data and covariates
        data, covar = self._prepare_data(X.copy())

        # get batch effect levels
        self.batch_levels = covar["SITE"].unique().tolist()

        # train ComBat model
        if not self.smooth_terms:

            self._vprint("++ No GAM terms specified; using linear model ++")

            if self.random_state: np.random.seed(self.random_state)

            self.model, _ = neuroHarmonize.harmonizationLearn(
                data.values,
                covar
            )

        else:
            
            self._vprint(f"++ Using GAM terms {self.smooth_terms} ++")

            if self.random_state: np.random.seed(self.random_state)

            self.model, _ = neuroHarmonize.harmonizationLearn(
                data.values,
                covar,
                smooth_terms = self.smooth_terms,
                smooth_term_bounds = (None, None) if (self.smooth_bounds is None) else self.smooth_bounds 
            )

        return self
    
    def transform(self, X, y = None):

        self._vprint("++ Applying ComBat to test data ++")

        # separate into data and covariates
        data, covar = self._prepare_data(X.copy())

        # apply ComBat model
        data_combat = neuroHarmonize.harmonizationApply(
            data.values,
            covar,
            self.model
        )

        # convert to pandas df and concatenate with all non-feature columns
        combat = pd.concat([
            X.copy().drop(self.features, axis = 1).reset_index(drop = True),
            pd.DataFrame(data_combat, columns = data.columns)
        ], axis = 1)

        return combat

    def get_param_dict(self):

        """
        Store GAM-ComBat model parameters in dictionary

        Parameters
        ----------
        None

        Returns
        -------
        dict
            dictionary of model parameters
        """
        
        # store params in dictionary
        param_dict = {
            "features": self.features.tolist(),
            "covariates": self.covariates,
            "site_name": self.site_name,
            "smooth_terms": self.smooth_terms,
            "smooth_bounds": self.smooth_bounds,
            "random_state": self.random_state
        }
        
        if hasattr(self, "batch_levels"):
            param_dict["batch_levels"] = self.batch_levels

        return param_dict


class ComBatCV():

    """
    Perform a repeated stratified (grouped) cross-validation to harmonize
    a dataset using GAM-ComBat. The procedure is roughly the following:

        1. Get stratified (grouped) k-fold split
        2. Choose 1 fold for hold-out test set, use remaining folds to
        train ComBat harmonization model
        3. Apply trained model to test set to harmonize
        4. Repeat steps 2 and 3 for remaining folds (inner loop) and
        concatenate harmonized test folds to get the full dataset
        5. Repeat steps 1 through 4 for number of repeats (outer loop)
        6. Average across repeats for the final harmonized dataset
    """

    def __init__(
        self,
        data_df,
        model,
        stratify_col = None,
        group_col = None,
        n_splits = 5,
        n_repeats = 20,
        random_state = None,
        verbose = False
    ):
        """
        Parameters
        ----------
        data_df : pandas.DataFrame, shape (n,m)
            dataframe containing data to harmonize; this includes features
            as well as covariates
        model : GAMComBat
            GAM-ComBat model instance; specify the features, covariates and
            smooth terms through this class rather than this function
        stratify_col : numpy.ndarray, shape (n,)
            column containing labels to use for stratifying CV splits
        group_col : numpy.ndarray, shape (n,)
            column containing labels to use for grouping rows in grouped CV;
            only used in method `cross_validate_grouped`
        n_splits : int
            number of CV splits
        n_repeats : int
            number of times to repeat CV
        random_state : int
            random state passed into sklearn.model_selection.RepeatedStratifiedKFold
            to obtain CV splits
        """

        self.data_df = data_df
        self.model = model
        self.stratify_col = stratify_col
        self.group_col = group_col
        self.n_splits = n_splits
        self.n_repeats = n_repeats
        self.random_state = random_state
        self.verbose = verbose

        # get batch levels
        _, _covar = self.model._prepare_data(self.data_df.copy())
        self.batch_levels = _covar["SITE"].unique().tolist()

    def _vprint(self, s):

        if self.verbose: print(s)

    def cross_validate(self):

        """
        Main method to perform repeated stratified CV to harmonize data using
        GAM-ComBat; see procedure in the class docstring

        Parameters
        ----------
        None

        Returns
        -------
        pandas.DataFrame, shape (n,m)
            dataframe containing harmonized features, as well as remaining
            columns (these remain intact)
        """

        # define K-fold split
        kfold_split = RepeatedStratifiedKFold(
            n_splits = self.n_splits,
            n_repeats = self.n_repeats,
            random_state = self.random_state
        ).split(X = self.data_df, y = self.stratify_col)

        # perform CV
        combat_cv_list = []
        for i, (train_idx, test_idx) in enumerate(kfold_split):

            current_rep = int(np.ceil((i+1) / self.n_splits))
            current_fold = (i % self.n_splits) + 1
            self._vprint(f"++ Repetition {current_rep}, CV split {current_fold} ++")

            cv_model_loop = clone(self.model)
            cv_model_loop.fit(X = self.data_df.iloc[train_idx, :]);
            combat_cv_loop = cv_model_loop.transform(self.data_df.iloc[test_idx, :])

            combat_cv_list.append(combat_cv_loop.loc[:, np.concatenate([self.model.features, ["idx"]])])

        # bind rows
        combat_cv = pd.concat(combat_cv_list, axis = 0)

        # group by idx, compute mean across idx
        combat_cv_mean = combat_cv.groupby("idx").mean().reset_index()

        # add missing columns
        combat_cv_mean = pd.merge(
            left = self.data_df.drop(self.model.features, axis = 1).reset_index(drop = True),
            right = combat_cv_mean,
            how = "outer",
            on = "idx"
        )

        return combat_cv_mean

    def cross_validate_grouped(self):

        """
        Perform repeated grouped stratified CV to harmonize data. This is
        useful for separating subjects between train and test splits, such
        that no subject appears in both splits at the same time, e.g. when
        your dataset has multiple measurements of the same subject at
        different time points

        Parameters
        ----------
        None
            
        Returns
        -------
        pandas.DataFrame, shape (n,m)
            dataframe containing harmonized features, as well as remaining
            columns (these remain intact)
        """

        # set random seed
        if self.random_state: np.random.seed(self.random_state)

        # start repeated loop
        combat_cv_list = []
        for i in range(self.n_repeats):

            # define K-fold split
            kfold_split = StratifiedGroupKFold(
                n_splits = self.n_splits,
                shuffle = True
            ).split(
                X = self.data_df,
                y = self.stratify_col,
                groups = self.group_col
            )

            # perform CV
            for j, (train_idx, test_idx) in enumerate(kfold_split):

                self._vprint(f"++ Repetition {i}, CV split {j} ++")

                cv_model_loop = clone(self.model)
                cv_model_loop.fit(X = self.data_df.iloc[train_idx, :]);
                combat_cv_loop = cv_model_loop.transform(self.data_df.iloc[test_idx, :])

                combat_cv_list.append(combat_cv_loop.loc[:, np.concatenate([self.model.features, ["idx"]])])

        # bind rows
        combat_cv = pd.concat(combat_cv_list, axis = 0)

        # group by idx, compute mean across idx
        combat_cv_mean = combat_cv.groupby("idx").mean().reset_index()

        # add missing columns
        combat_cv_mean = pd.merge(
            left = self.data_df.drop(self.model.features, axis = 1).reset_index(drop = True),
            right = combat_cv_mean,
            how = "outer",
            on = "idx"
        )

        return combat_cv_mean

    def get_param_dict(self):

        """
        Store GAM-ComBat model parameters + CV parameters
        in dictionary

        Parameters
        ----------
        None

        Returns
        -------
        dict
            dictionary of model parameters
        """

        param_dict = self.model.get_param_dict()
        param_dict["batch_levels"] = self.batch_levels

        cv_param_dict = {
            "stratify_col_levels": pd.unique(self.stratify_col).tolist(),
            "grouped": False if self.group_col is None else True,
            "n_splits": self.n_splits,
            "n_repeats": self.n_repeats,
            "random_state": self.random_state
        }
        param_dict["cv_params"] = cv_param_dict

        return param_dict


class PETDataLoader:
    
    def __init__(self, data_dir):

        self.data_dir = data_dir

    @staticmethod
    def _read_rds(rds_path):

        # https://stackoverflow.com/questions/40996175/loading-a-rds-file-in-pandas

        from pyreadr import read_r

        df = read_r(rds_path)[None]

        return df

    @staticmethod
    def _create_dummy_variables(df, cols):

        dummy_df = pd.get_dummies(
            df.loc[:, cols],
            columns = cols
        )
        df = pd.concat([df, dummy_df], axis = 1)

        return df

    def load_data(self, all_scans = False, pvc = False, batch_other = None, white_matter = False):

        """
        Load PET data for use in ComBat harmonization and
        predictive modeling; default is to load only baseline
        scans but there's also the option of loading all
        PET scan data

        Parameters
        ----------
        all_scans : bool
            If True, load all PET scan data
        pvc : bool
            If True, only select scans that have been partial
            volume corrected; note that this will only select
            OASIS scans, since ADNI does not do PVC
        """

        def mutate_df(df):

            # add batch effect column named SITE
            if not batch_other is None:
                # add an additional batch effect
                df["SITE"] = df["study"] + "_" + df["tracer"] + "_" + df[batch_other].astype(str)
            else:
                # just use study and tracer
                df["SITE"] = df["study"] + "_" + df["tracer"]

            # create dummy variables for sex, apoe and clinical group
            cat_col = ["sex", "apoe", "clinical_group", "clinical_group_extended"]
            df = self._create_dummy_variables(df, cat_col)

            return df

        def get_train_df(data_df, crossover_df, pvc):

            # select all OASIS scans that are not part of the crossover set and use for training
            train_df = data_df[
                np.logical_and.reduce([
                    data_df["study"] == "OASIS",
                    data_df["pvc"] == pvc,
                    ~data_df["idx"].isin(np.concatenate([crossover_df["idx.AV45"], crossover_df["idx.PIB"]]))
                ])
            ]

            return train_df

        wm_suffix = "_wm" if white_matter else ""

        # load data
        if pvc:
            self.crossover = self._read_rds(path.join(self.data_dir, "crossover_pvc.RDS"))
            if all_scans:
                self.data = self._read_rds(path.join(self.data_dir, f"pet{wm_suffix}.RDS"))
                self.data = self.data[self.data["pvc"]]
            else:
                self.data = self._read_rds(path.join(self.data_dir, "baseline_pvc.RDS"))
        else:
            self.crossover = self._read_rds(path.join(self.data_dir, f"crossover{wm_suffix}.RDS"))
            if all_scans:
                self.data = self._read_rds(path.join(self.data_dir, f"pet{wm_suffix}.RDS"))
                self.data = self.data[~self.data["pvc"]]
            else:
                self.data = self._read_rds(path.join(self.data_dir, f"baseline{wm_suffix}.RDS"))

        self.feature_names = np.loadtxt(path.join(self.data_dir, "csv", "feature_names_trunc.csv"), dtype = str)
        self.wm_feature_names = np.loadtxt(path.join(self.data_dir, "csv", "white_matter_regions.csv"), dtype = str)

        # add new columns
        self.data = mutate_df(self.data)
        self.crossover = mutate_df(self.crossover)

        # drop nan rows (missing covariates)
        self.data = self.data.dropna(subset = ["age", "sex", "apoe"])
        
        # get training data for crossover experiment
        self.crossover_train = get_train_df(self.data, self.crossover, pvc)

        # separate amyloid +/-
        self.amypos = self.data[self.data["amyloid_positive"]]
        self.amyneg = self.data[~self.data["amyloid_positive"]]

        return None


# ======== helper functions ========

def get_smooth_bounds(df_list, col):

    v = np.concatenate([df[col] for df in df_list])
    return (v.min(), v.max())

def dict2json(d, filepath):

    with open(filepath, "w") as fp:
        json.dump(d, fp, indent = 4)

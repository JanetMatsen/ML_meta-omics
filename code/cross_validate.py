import os
import pandas as pd
import re

import matplotlib.pyplot as plt

from CCA import ExpressionCCA

class CrossValCCA(object):
    def __init__(self, raw_data_path, uv_dir,
                 input_filepath,
                 path_to_R_script='../../code/sparse_CCA.R',
                 ):
        self.raw_data_path = raw_data_path
        self.uv_dir = uv_dir
        self.input_filepath = input_filepath
        self.path_to_R_script=path_to_R_script

        self.results = pd.DataFrame()
        self.models = dict()
        self.models_made = 0

    def model(self, x_train_filepath, z_train_filepath, pen_x, pen_z,
              x_val_filepath, z_val_filepath,
              x_gene_filepath, z_gene_filepath,  # gene names
              verbose=False):

        self.models_made += 1
        # todo: remove expected filename if it exsits.
        cca = ExpressionCCA(x_train_filepath= x_train_filepath,
                            z_train_filepath= z_train_filepath,
                            penalty_x = pen_x,
                            penalty_z = pen_z,
                            x_val_filepath= x_val_filepath,
                            z_val_filepath= z_val_filepath,
                            x_gene_filepath=x_gene_filepath,
                            z_gene_filepath=z_gene_filepath,
                            input_filepath = self.input_filepath,
                            u_v_output_dir = self.uv_dir,
                            verbose = verbose,
                            path_to_R_script=self.path_to_R_script)
        # Store the model in the object's model dict
        self.models[self.models_made] = cca

        # build up the objects results df
        cca.summarise()
        model_summary = cca.summary
        model_summary['model number'] = self.models_made
        model_summary['penalty_x'] = pen_x
        model_summary['penalty_z'] = pen_x
        model_summary['fold'] = self.fold_name(x_train_filepath)
        model_summary['model'] = cca

        model_row = pd.DataFrame(self.prep_for_pandas(model_summary))
        self.results = pd.concat([self.results, model_row], axis = 0)

    @staticmethod
    def prep_for_pandas(mydict):
        return pd.DataFrame({k:[v] for k, v in mydict.items()})

    @staticmethod
    def fold_name(string):
        m = re.search('[_A-z]+fold([0-9]+)[._A-z]+', string)
        return int(m.group(1))

    def plot_series_by_fold(self, x_col, y_col, title,
                            colors=None, figsize=(4, 4), filename=None):
        fig, ax = plt.subplots(1, 1, figsize=figsize)

        if colors is 'greens':
            colors = ['#74c476', '#41ab5d', '#238b45', '#005a32'] # http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=5
        elif colors is 'purples':
            colors = ['#8c96c6', '#8c6bb1', '#88419d', '#6e016b']

        for tup, df in self.results.groupby('fold'):
            fold = int(tup)
            df = df.sort(x_col)
            plt.plot(df[x_col], df[y_col], color=colors[fold-1],
                     linestyle='--', marker='o', label="fold {}".format(fold))
        plt.legend(loc = 'best')
        plt.xlabel(x_col)
        plt.ylabel(y_col)
        if title is not None:
            plt.title(title)

        plt.tight_layout()

        if filename is not None:
            fig.savefig(filename)

        return fig

    def plot_correlation_vs_penalty(self, set, penalty='x',
                                    figsize=(4, 4), title=None, filename=None):

        if penalty == 'x':
            x = 'penalty_x'
        elif penalty == 'z':
            x = 'penalty_z'
        else:
            raise NameError, "did you mean 'x' or 'z'?"

        if set == 'train':
            y = 'train correlation'
            title = "cross-validation fit: training data"
        elif set == 'val' or set =='validation':
            y = 'validation correlation'
            title = "cross-validation fit: validation data"
        else:
            raise NameError, "did you mean 'train' or 'val'?"

        fig = self.plot_series_by_fold(x_col=x, y_col=y, title=title,
                                        colors='purples', figsize=figsize,
                                        filename=filename)
        return fig

    def plot_num_nonzero_coeffs_vs_penalty(self, set, title=None,
                                           figsize=(4, 4), filename=None):

        if set == 'x':
            x = 'penalty_x'
            y = '# nonzero u weights'
        if set == 'z':
            x = 'penalty_z'
            y = '# nonzero v weights'

        fig = self.plot_series_by_fold(x_col=x, y_col=y, title=title,
                                        colors='purples', figsize=figsize,
                                        filename=filename)

        return fig





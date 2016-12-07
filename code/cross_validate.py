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
              x_val_filepath, z_val_filepath, gene_name_filepath,
              verbose=False
              ):

        self.models_made += 1
        # todo: remove expected filename if it exsits.
        cca = ExpressionCCA(x_train_filepath= x_train_filepath,
                            z_train_filepath= z_train_filepath,
                            penalty_x = pen_x,
                            penalty_z = pen_z,
                            x_val_filepath= x_val_filepath,
                            z_val_filepath= z_val_filepath,
                            gene_filepath = gene_name_filepath,
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

    def plot_sparsity_projections(self, set='x', figsize=(4, 4), filename=None):
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        if set == 'x':
            x = 'penalty_x'
            y = '# nonzero u weights'
        if set == 'z':
            x = 'penalty_z'
            y = '# nonzero v weights'
        for tup, df in self.results.groupby('fold'):
            fold = int(tup)
            print(fold)
            colors = ['#8c96c6', '#8c6bb1', '#88419d', '#6e016b']
            df = df.sort(x)
            plt.plot(df[x], df[y], color=colors[fold-1],
                     linestyle='--', marker='o', label="fold {}".format(fold))
        plt.legend(loc = 'best')
        plt.xlabel(x)
        plt.ylabel('# nonzero weights')

        if filename is not None:
            fig.savefig(filename + '.pdf')

        return fig





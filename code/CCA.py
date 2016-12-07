import matplotlib.pylab as plt
import numpy as np
import os
import re
import subprocess

import pandas as pd

class CcaAnalysis(object):
    """
    Analyze the results of sparse CCA results (run in R)
    """
    def __init__(self, x, z, u, v, val_x, val_z):
        self.x = x
        self.z = z
        self.u = u
        self.v = v
        self.val_x = val_x
        self.val_z = val_z

        self.project()
        self.summary = None

    def project(self):
        self.x_projected = self.x.dot(self.u)
        self.x_val_projected = self.val_x.dot(self.u)
        self.z_projected = self.z.dot(self.v)
        self.z_val_projected = self.val_z.dot(self.v)

    def plot_projections(self, filename=None):
        # scatter plot of x.dot(u) vs z.dot(v)
        # one color for train, one color for validation
        fig, ax = plt.subplots(1, 1, figsize=(3,3))
        colors = ['#bdbdbd','#31a354']
        plot_vars = [(self.x_projected, self.z_projected),
                    (self.x_val_projected, self.z_val_projected)]

        series = 0
        for (x, y) in plot_vars:
            plt.scatter(x, y, linestyle='--', marker='o', color=colors[series])
            series += 1
        plt.legend(loc = 'best')

        if filename is not None:
            fig.savefig(filename + '.pdf')

#        print('testing')
        return fig

    @staticmethod
    def correlation(x, z):
        # return correlation of projections for x.dot(u) and z.dot(v)
        corr_matrix = np.corrcoef(x,z)
        assert corr_matrix.shape == (2,2)
        return corr_matrix[0,1]

    def summarise(self):
        summary = {}
        summary['train correlation'] = \
            self.correlation(self.x_projected, self.z_projected)
        summary['validation correlation'] = \
            self.correlation(self.x_val_projected, self.z_val_projected)
        summary['# nonzero u weights'] = self.num_nonzero(self.u)
        summary['# nonzero v weights'] = self.num_nonzero(self.v)

        #summary = {k:[v] for k, v in summary.items()}
        #return pd.DataFrame(summary)
        self.summary = summary
        
    def get_summary(self):
        if self.summary is None:
            self.summarise()
        return self.summary
            

    @staticmethod
    def num_nonzero(vector):
        # for counting # of nonzero components in u and v
        return sum(vector != 0)


class ExpressionCCA(CcaAnalysis):
    """
    Wrapper class: prepare and run CCA for particular x, z data sets
    """
    def __init__(self, x_train_filepath, z_train_filepath,
                 x_val_filepath, z_val_filepath,
                 gene_filepath,
                 input_filepath, u_v_output_dir,
                 penalty_x, penalty_z,
                 verbose = False,
                 path_to_R_script='../../code/sparse_CCA.R'):

        self.penalty_x = penalty_x
        self.penalty_z = penalty_z
        self.path_to_R_script=path_to_R_script

        assert os.path.exists(input_filepath), \
            "input filepath, {}, doesn't exist".format(input_filepath)
        if not os.path.exists(u_v_output_dir):
            os.mkdir(u_v_output_dir)
        self.u_v_output_dir = u_v_output_dir

        self.x_train_filepath = x_train_filepath
        self.z_train_filepath = z_train_filepath
        self.x_val_filepath = x_val_filepath
        self.z_val_filepath = z_val_filepath
        self.gene_filepath = gene_filepath

        self.penalty_x = penalty_x
        self.penalty_z = penalty_z

        # prepare u and v
        x, z, u, v = self.write_csv_and_run_R(verbose=verbose)

        super(ExpressionCCA, self).__init__(
            x=x, z=z, u=u, v=v,
            val_x=self.load_array(self.x_val_filepath),
            val_z=self.load_array(self.z_val_filepath))

    @staticmethod
    def load_array(filepath, verbose=False):
        vector = np.genfromtxt(filepath, delimiter='\t')
        if verbose:
            print('vector {} has shape {}'.format(filepath, vector.shape))
        return vector

    def write_csv_and_run_R(self, delete_u_v=False, verbose=False):
        x = self.load_array(self.x_train_filepath)
        self.x = x
        z = self.load_array(self.z_train_filepath)
        self.z = z
        names = self.load_array(self.gene_filepath)
        self.gene_names = names

        # get the data back out
        def prepare_output_filename(input_filename, extra_string):
            # methylotroph_fold1_train.tsv --> fold1_train_u.tsv
            if verbose:
                print('output dir: {}'.format(self.u_v_output_dir))
            s = os.path.basename(input_filename)
            m = re.search('[_A-z]+(fold[0-9]+[._A-z]+.tsv)', s)
            s = m.group(1)
            s = s.replace('.tsv', '_{}.tsv'.format(extra_string))
            s = os.path.join(self.u_v_output_dir, s)
            if verbose:
                print('Will save output for {} to {}'.format(input_filename, s))
            return s

        u_path = prepare_output_filename(self.x_train_filepath,
                                         extra_string='u_penX' + str(self.penalty_x)
                                        + '_penZ' + str(self.penalty_z))
        v_path = prepare_output_filename(self.z_train_filepath ,
                                         extra_string='v_penX' + str(self.penalty_x)
                                        + '_penZ' + str(self.penalty_z))

        # Run R
        # todo: this will keep concatenating to the same file.  Need to delete
        # it sometimes, put a time stamp on it, or something else.
        stdout_file = open('stdout_CCA.txt', 'a')
        stderr_file = open('stderr_CCA.txt', 'a')
        command = ['Rscript', self.path_to_R_script,
                   self.x_train_filepath, self.z_train_filepath,
                   u_path, v_path,
                   str(self.penalty_x), str(self.penalty_z)]
        if verbose:
            print('command: \n {}'.format(" ".join(command)))
        subprocess.check_call(command, stdout=stdout_file, stderr=stderr_file)
        stdout_file.close()
        stderr_file.close()

        # R adds a header row, 'V1' we chop off.
        print('read u in from {}'.format(u_path))
        print('read v in from {}'.format(v_path))
        u = np.genfromtxt(u_path, delimiter='\t', skip_header=1)
        v = np.genfromtxt(v_path, delimiter='\t', skip_header=1)
        if verbose:
            print(u.shape)
            print(v.shape)

        # todo: assert some shape constraints.

        return x, z, u, v


class SklearnCca(CcaAnalysis):
    def __init__(self, penalty_x, penalty_z,
                 output_dir = './sklearn_test',
                 path_to_R_script='../../code/sparse_CCA.R'):
        # Recreate the CCA results from sklearn
        # http://scikit-learn.org/stable/auto_examples/cross_decomposition/plot_compare_cross_decomposition.html#sphx-glr-auto-examples-cross-decomposition-plot-compare-cross-decomposition-py

        self.penalty_x = penalty_x
        self.penalty_z = penalty_z
        self.path_to_R_script=path_to_R_script
        self.output_dir = output_dir
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        X_train, Y_train, X_test, Y_test = self.prepare_data()

        # populate the
        u, v = self.write_csv_and_run_R(X_train, Y_train)

        # Fill in the CCA analysis class instance
        super(SklearnCca, self).__init__(x=X_train, z=Y_train,
                                         u=u, v=v,
                                         val_x=X_test, val_z=Y_test)

    def prepare_data(self):
        """
        Pasted directly from the sklearn example:
        http://scikit-learn.org/stable/auto_examples/cross_decomposition/plot_compare_cross_decomposition.html#sphx-glr-auto-examples-cross-decomposition-plot-compare-cross-decomposition-py
        Note: their X, Y in Sklearn is x, z in R's sparse CCA
        """
        n = 500
        # 2 latents vars:
        l1 = np.random.normal(size=n)
        l2 = np.random.normal(size=n)

        latents = np.array([l1, l1, l2, l2]).T
        # X, Y are 250 examples with 4 features each.
        X = latents + np.random.normal(size=4 * n).reshape((n, 4))
        Y = latents + np.random.normal(size=4 * n).reshape((n, 4))

        X_train = X[:n / 2]
        Y_train = Y[:n / 2]
        X_test = X[n / 2:]
        Y_test = Y[n / 2:]

        # Python has column features and row samples
        # R has row samples and column features.
        #return X_train.T, Y_train.T, X_test.T, Y_test.T
        return X_train, Y_train, X_test, Y_test

    def write_csv_and_run_R(self, X_train, Y_train, delete_u_v=False):
        print(X_train.shape)
        print(Y_train.shape)
        arrays = [X_train, Y_train]
        names = ['X', 'Y']
        filenames_to_remove = []

        for a, n in zip(arrays, names):
            fname = n + '.tsv'
            # TODO: remove filename if it exists (precaution for script failure)
            filenames_to_remove.append(fname)
            np.savetxt(fname, a, delimiter='\t')

        u_path = os.path.join(self.output_dir, 'u.tsv')
        v_path = os.path.join(self.output_dir, 'v.tsv')

        # Run R
        command = ['Rscript', self.path_to_R_script,
                   filenames_to_remove[0], filenames_to_remove[1],
                   u_path, v_path,
                   str(self.penalty_x), str(self.penalty_z)]
        print('command: \n {}'.format(" ".join(command)))

        file_handle = open('stdout_sklearn.txt', 'a')
        subprocess.check_call(command, stdout=file_handle)

        u = np.genfromtxt(u_path, delimiter='\t', skip_header=1)
        v = np.genfromtxt(v_path, delimiter='\t', skip_header=1)
        print(u.shape)
        print(v.shape)

        # delete the files
        if delete_u_v:
            for filename in filenames_to_remove:
                if os.path.exists(filename):
                    os.remove(filename)

        return u, v

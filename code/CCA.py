import matplotlib.pylab as plt
import numpy as np
import os
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

        print('testing')
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

    @staticmethod
    def num_nonzero(vector):
        # for counting # of nonzero components in u and v
        return sum(vector == 0)



class SklearnCca(CcaAnalysis):
    def __init__(self, penalty_x, penalty_z,
                 path_to_R_script='../../code/sparse_CCA.R'):
        # Recreate the CCA results from sklearn
        # http://scikit-learn.org/stable/auto_examples/cross_decomposition/plot_compare_cross_decomposition.html#sphx-glr-auto-examples-cross-decomposition-plot-compare-cross-decomposition-py

        self.penalty_x = penalty_x
        self.penalty_z = penalty_z
        self.path_to_R_script=path_to_R_script

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

        # Run R
        command = ['Rscript', self.path_to_R_script,
                   filenames_to_remove[0], filenames_to_remove[1],
                   str(self.penalty_x), str(self.penalty_z), '.']
        print('command: \n {}'.format(" ".join(command)))
        subprocess.check_call(command)

        # get the data back out
        u_path = filenames_to_remove[0].replace('.tsv', '_u.tsv')
        v_path = filenames_to_remove[1].replace('.tsv', '_v.tsv')
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

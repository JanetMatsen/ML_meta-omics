import seaborn as sns
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA


class GenePCA():
    def __init__(self, raw_x, sample_names, gene_names, sample_info):
        """

        :param raw_x:
        :param colnames:
        :param sample_info:
        """
        self.pca = PCA()
        self.X = raw_x
        self.sample_names = sample_names
        self.gene_names = gene_names
        self.sample_info = sample_info

    def fit(self):
        # do PCA
        self.pca.fit(self.X)
        self.explained_variance = self.pca.explained_variance_ratio_
        self.X_transformed = self.pca.transform(self.X)


        # prepare data for plotting
        plot_data = pd.DataFrame({'direction 1': self.X_transformed[:, 0],
                                  'direction 2': self.X_transformed[:, 1]})
        # append on sample names
        plot_data = pd.concat([plot_data, self.sample_names], axis=1)
        # append on the data descriptors
        plot_data = plot_data.merge(self.sample_info,
                                    left_on = 'sample', right_on = 'ID')
        self.plot_data = plot_data

        # prepare a dataframe with the most differentiating genes.
        component_features = pd.DataFrame(self.pca.components_.T)
        component_features['abs sum'] = component_features.abs().sum(axis=1)
        component_features = pd.concat([component_features,
                                        self.gene_names], axis=1)
        component_features.sort(columns='abs sum', ascending=False,
                                inplace=True)
        self.component_features = component_features

        print(self.plot_data.head())


def plot_pca_results(plot_data, variances,
                     color_palette,
                     color_variable,
                     facet_row=True, uniform_axes=True,
                     main_dir='./', plot_dir='./figures/',
                     savefig=False, subplot_size=3,
                     filename=None):

    # prepare axis labels, which also serve as dataframe column names.
    x_axis_label = 'principal component 1 ({0:0.2%})'.format(variances[0])
    y_axis_label = 'principal component 2 ({0:0.2%})'.format(variances[1])
    plot_data = plot_data.rename(columns={'direction 1':x_axis_label})
    plot_data = plot_data.rename(columns={'direction 2':y_axis_label})


    # update matplotlib params for bigger fonts, ticks:
    mpl.rcParams.update({
        'font.size': 16, 'axes.titlesize': 17, 'axes.labelsize': 15,
        'xtick.labelsize': 10, 'ytick.labelsize': 13,
        'font.weight': 600,
        'axes.labelweight': 600, 'axes.titleweight': 600})

    # Plot with Seaborn
    if facet_row:
        plt.figure(figsize=(4, 8))
    else:
        plt.figure(figsize=(6, 12))
    sns.set(style="ticks")

    # prepare the max and min axes values if we are forcing them to same range
    pc_colnames = [col for col in plot_data.columns
                   if 'principal component' in col]

    max_value = plot_data[pc_colnames].max(axis=0).max()
    min_value = plot_data[pc_colnames].min(axis=0).min()

    axis_max = math.ceil(max_value * 100) / 100.0
    axis_min = math.floor(min_value * 100) / 100.0

    def base_plot(**kwargs):
        plot = sns.FacetGrid(plot_data,
                             hue=color_variable, palette=color_palette,
                             size=subplot_size, aspect=1,
                             **kwargs)
        plot = (plot.map(plt.scatter, x_axis_label, y_axis_label,
                         edgecolor="w", s=60).add_legend())
        return plot

    plot_args = {}

    if facet_row:
        plot_args['row'] = 'oxy'
        plot_args['col'] = 'rep'
    if uniform_axes:
        plot_args['xlim'] = (axis_min, axis_max)
        plot_args['ylim'] = (axis_min, axis_max)

    if len(plot_args) > 0:
        print(plot_args)

    g = base_plot(**plot_args)

    if filename is None and savefig:
        filename = 'plot'



    if savefig:
        # prepare a filename, depending on whether all taxonomy or only genus
        # is used.
        if uniform_axes:
            filename += '_unif_axes_'
        if facet_row:
            filename += '--faceted.pdf'
        else:
            filename += '.pdf'
        g.fig.savefig(filename)
    else:
        return g

def plot_color_by_week(plot_data, variances, plot_pca_kwargs):
    # define a custom color palette using:
    # Conditions were switched at week ten, so seven early samples in
    # the original condition and four latest samples in an alternative
    # condition.
    color_palette = build_double_color_palette(num_items=14 - 4 + 1,
                                               weeks_before_switch=7)

    return plot_pca_results(plot_data, variances,
                            color_palette=color_palette,
                            color_variable='week',
                            **plot_pca_kwargs)


def plot_color_by_category(plot_data, variances,
                           categories, n_categories,
                           plot_pca_kwargs):
    color_palette = build_categorical_color_palette(n_categories, )
    return plot_pca_results(plot_data, variances,
                            color_palette=color_palette,
                            color_variable=categories,
                            **plot_pca_kwargs)


def build_categorical_color_palette(n):
    return sns.color_palette("Set1", n)


def build_double_color_palette(num_items, weeks_before_switch):
    num_pre_switch_colors = weeks_before_switch
    num_post_switch_colors = num_items - num_pre_switch_colors
    print('preparing colors for {} pre-oxygen-switch'.format(
        num_pre_switch_colors),
          'samples and {} post-switch samples'
          .format(num_post_switch_colors))

    # get the first colors from this pallete:
    pre_switch_colors = \
        sns.cubehelix_palette(11, start=.5, rot=-.75)[0:num_pre_switch_colors]
    print(pre_switch_colors)

    # get post-switch colors here:
    # post_switch_colors = sns.diverging_palette(220, 20,
    # n=6)[::-1][0:num_post_switch_colors]
    post_switch_colors = \
        sns.color_palette("coolwarm", num_post_switch_colors)
    # sns.light_palette("navy", reverse=True)[0:num_post_switch_colors]
    rgb_colors = pre_switch_colors + post_switch_colors
    sns.palplot(rgb_colors)

    # check that we got the right amount
    print(num_items)
    assert (num_items == len(rgb_colors))
    print("")
    return rgb_colors
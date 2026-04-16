import glob
import math
import os
import unittest
from enum import auto
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt, image as mpimg
import config
from figure_generation import ks_modeling
from figure_generation.coalescent_plot_aggregation import get_run_time_in_minutes, read_data_csv
from figure_generation.histogram_plotter import read_Ks_csv
from figure_generation.multi_row_Ks_aggregation_plots import SPKSvsDGKSPlotData


class AUTOvsALLOPlotData:
    def __init__(self, auto_out_path, auto_run_list, auto_run_list_name,
                 allo_out_path, allo_run_list,
                 bin_sizes_Ks,
                 xmax_Ks, ymax_Ks,
                 show_KS_predictions,
                 include_annotation, which_plot_panels_to_show_legend, plot_title_lamda):
        self.auto_out_path = auto_out_path
        self.auto_run_list = auto_run_list
        self.auto_run_list_name =auto_run_list_name
        self.allo_out_path = allo_out_path
        self.allo_run_list = allo_run_list
        self.bin_sizes_Ks = bin_sizes_Ks
        self.xmax_Ks = xmax_Ks
        self.ymax_Ks = ymax_Ks
        self.show_KS_predictions = show_KS_predictions
        self.include_annotation = include_annotation
        self.which_plot_panels_to_show_legend = which_plot_panels_to_show_legend
        self.plot_title_lamda = plot_title_lamda

def plot_allo_vs_auto_ks(this_ax, allo_config_used, allo_ks_by_gene,
                         auto_config_used, auto_ks_by_gene,
                         title, bin_size, xmax, ymax, show_predictions,
                         include_annotation, show_legend):

    num_allo_genes = len(allo_ks_by_gene)
    num_auto_genes = len(auto_ks_by_gene)

    if not xmax:
        xmax = max(allo_ks_by_gene)
    bins = np.arange(0, xmax, bin_size)
    dgx_hist_ys=[]

    if len(allo_ks_by_gene) > 0:

        my_label = 'Allo Ks by gene'
        if include_annotation:
            my_label = 'Allo Ks by gene\n'\
                           + "(" + str(num_allo_genes) + " paralogs in genome)"
        dgx_hist_ys, bins, patches = this_ax.hist(allo_ks_by_gene, bins=bins, facecolor='b', alpha=0.25,
                                                  label=my_label,
                                                  density=False)

    if len(auto_ks_by_gene) > 0:

        my_label = 'Auto Ks by gene'
        if include_annotation:
            my_label = 'Auto Ks by gene\n'\
                           + "(" + str(num_auto_genes) + " paralogs in genome)"
        this_ax.hist(auto_ks_by_gene, bins=bins, facecolor='g', alpha=0.50,
                     label=my_label,
                     density=False)

    #this_ax.axvline(x=allo_config_used.t_div_as_ks, color='b', linestyle='--', label="input Tdiv as Ks")
    if allo_config_used:
        theoretical_ks_mean_now= allo_config_used.mean_Ks_from_Tc + allo_config_used.t_div_as_ks
        this_ax.axvline(x=theoretical_ks_mean_now, color='purple', linestyle='--',
                    label="Expected Ks mean (Allo)",lw=3)
    if auto_config_used:
        this_ax.axvline(x=auto_config_used.mean_Ks_from_Nb, color='k', linestyle='--',
                    label="Expected Ks mean (Auto)",lw=3)

    #mean_Ks_from_Nb_string=  "({:.2E})".format(config_used.mean_Ks_from_Nb)
    #this_ax.axvline(x=config_used.mean_Ks_from_Nb,
    #                color='g', linestyle='--', label="Tc due to Nb " + mean_Ks_from_Nb_string )

    if allo_config_used:
        add_allo_auto_Ks_annotations(this_ax, allo_config_used, allo_ks_by_gene, dgx_hist_ys, bins,
                                 include_annotation, show_predictions)
    else:
        add_allo_auto_Ks_annotations(this_ax, auto_config_used, auto_ks_by_gene, dgx_hist_ys, bins,
                                 include_annotation, show_predictions)
    if ymax:
        this_ax.set(ylim=[0, ymax])
    this_ax.set(xlim=[0, xmax])
    this_ax.set(xlabel="Ks")
    this_ax.set(title=title)

    if show_legend:
        this_ax.legend()


def make_multirow_Allo_vs_Auto_fig_with_subplots(mulitplotdatalist, png_out):

    num_plot_rows=len(mulitplotdatalist)
    num_runs = len(mulitplotdatalist[0].auto_run_list)
    par_dir = Path(__file__).parent.parent
    image_folder = os.path.join(par_dir, "images")
    png_Tnow = os.path.join(image_folder, 'Auto_vs_Allo_now_crop_2.png')

    if mulitplotdatalist[0].include_annotation:
        fig, ax = plt.subplots(num_plot_rows, num_runs, figsize=(40, 25))
        dpi_needed = 350
    else:
        fig, ax = plt.subplots(num_plot_rows, num_runs, figsize=(20,num_plot_rows*3))
        dpi_needed = 100

    captions=[
            "polyploid Ks at T_now\nfor Autopolyploid",
            "polyploid Ks at T_now\nfor Allo and Autopolyploid",
            "ancestral Tc at T_div\nfor Allo and T_wgd for Auto"]

    for plot_row_index in range(0,len(mulitplotdatalist)):

        mulitplotdata = mulitplotdatalist[plot_row_index]
        fig.suptitle(mulitplotdata.auto_run_list_name)

        col_offset_for_cartoon=0
        if not mulitplotdata.auto_run_list[0]:
            plot_expository_allo_auto_image(ax, png_Tnow, captions[1],plot_row_index)
            col_offset_for_cartoon = 1

        for i in range(col_offset_for_cartoon, len(mulitplotdata.allo_run_list)):

            allo_run_name = mulitplotdata.allo_run_list[i]
            if allo_run_name:

                allo_run_path = os.path.join(mulitplotdata.allo_out_path, allo_run_name)
                print("dgx_run_path: " +allo_run_path )
                glob_results=glob.glob(allo_run_path + '/*.used.xml')
                print(glob_results)
                allo_input_xml_file = glob_results[0]
                allo_config_used = config.DemographiKS_config(allo_input_xml_file)
                csv_file_name = allo_config_used.sim_name + '.csv'
                ks_file = os.path.join(allo_run_path, csv_file_name)
                print("reading " + ks_file)
                allo_ks_results = read_Ks_csv(ks_file, False)
                allo_run_duration_in_m, allo_version = get_run_time_in_minutes(allo_run_path)

                plot_title = str(mulitplotdata.plot_title_lamda(allo_config_used))
                if mulitplotdata.include_annotation:
                    plot_title = "Ks at Tnow\n" + "burnin time=" + str(allo_config_used.burnin_time) + " gen,\n" \
                         + "Na=" + str(allo_config_used.ancestral_Ne)  + ", Nb=" + str(allo_config_used.bottleneck_Ne) +\
                             ", Tdiv=" + str(allo_config_used.DIV_time_Ge) + ", RC=" + str(allo_config_used.recombination_rate)


            else:
                allo_config_used = False
                allo_ks_results = []
                allo_run_duration_in_m = 0
                plot_title = "Ks at Tnow"


            auto_run_name = mulitplotdata.auto_run_list[i]
            if auto_run_name:

                auto_run_path = os.path.join(mulitplotdata.auto_out_path,auto_run_name)
                print("auto_run_path: " +auto_run_path )
                print(auto_run_path)
                glob_results=glob.glob(auto_run_path + '/*.used.xml')
                auto_input_xml_file = glob_results[0]
                auto_config_used = config.DemographiKS_config(auto_input_xml_file)
                csv_file_name = auto_config_used.sim_name + '.csv'
                ks_file = os.path.join(auto_run_path, csv_file_name)
                print("reading " + ks_file)
                auto_ks_results = read_Ks_csv(ks_file, False)
                auto_run_duration_in_m, auto_version = get_run_time_in_minutes(auto_run_path)

                if not allo_run_name:
                    plot_title = "Ks at Tnow\n" + str(mulitplotdata.plot_title_lamda(auto_config_used))
            else:
                spx_ks_results = []
                spx_run_duration_in_m = 0
                specks_mrcas_by_gene = False

            show_legend = (i in mulitplotdata.which_plot_panels_to_show_legend)
            plot_allo_vs_auto_ks(ax[plot_row_index, i], allo_config_used, allo_ks_results,
                                 auto_config_used, auto_ks_results,
                    plot_title, mulitplotdata.bin_sizes_Ks[i],
                                 mulitplotdata.xmax_Ks[i],
                                 mulitplotdata.ymax_Ks[i],
                                 mulitplotdata.show_KS_predictions,
                                 mulitplotdata.include_annotation,show_legend)

            if allo_run_name:
                allo_Tc_csv_file = os.path.join(allo_run_path, "1_5_simulated_ancestral_gene_mrcas.csv")
                if not os.path.exists(allo_Tc_csv_file):
                    allo_Tc_csv_file = os.path.join(allo_run_path, "simulated_ancestral_gene_mrcas.csv")

                loci, allo_mrcas_by_gene = read_data_csv(allo_Tc_csv_file)

                plot_title = "Tcoal at Tdiv\nburnin time=" + str(auto_config_used.burnin_time) + " gen, " \
                         + "Na=" + str(auto_config_used.ancestral_Ne)

            auto_mrcas_by_gene=[]
            if auto_run_name:
                auto_Tc_csv_file = os.path.join(auto_run_path, "1_5_simulated_ancestral_gene_mrcas.csv")
                if not os.path.exists(auto_Tc_csv_file):
                    auto_Tc_csv_file = os.path.join(auto_run_path, "simulated_ancestral_gene_mrcas.csv")

                loci, auto_mrcas_by_gene = read_data_csv(auto_Tc_csv_file)

                plot_title = "Tcoal at Tdiv\nburnin time=" + str(auto_config_used.burnin_time) + " gen, " \
                         + "Na=" + str(auto_config_used.ancestral_Ne)

            theory_mrcas_by_gene=False
            if allo_config_used:
                ancestral_Ne = allo_config_used.ancestral_Ne
                num_genes =  allo_config_used.num_genes
            else:
                ancestral_Ne = auto_config_used.ancestral_Ne
                num_genes =  auto_config_used.num_genes


    ax[0, 1].set(ylabel="# paralog pairs in bin")
    ax[1, 1].set(ylabel="# genes in bin")
    plt.tight_layout()
    plt.savefig(png_out, dpi=dpi_needed)
    plt.clf()
    plt.close()

def add_allo_auto_Ks_annotations(this_ax, config_used, dgks_Ks_results, dgks_hist_results,
                                 bins, include_annotation, show_predictions):

    ks_predictions = ks_modeling.Ks_modeling_predictions(config_used, bins)
    if len(dgks_Ks_results) <= 0:
        return

    ks_fits= ks_modeling.Ks_modeling_fits(dgks_Ks_results, dgks_hist_results, bins)

    #compare simulated and predicted Ks means
    theoretical_ks_mean_now_as_string="Expected Ks mean ({:.2E})".format(ks_predictions.theoretical_ks_mean_now)
    simulated_ks_mean_now_as_string="Simulated Ks mean ({:.2E})".format(ks_fits.mean_ks_now_from_slim)

    #compare simulated and predicted Ks sigma
    theoretical_sigma_from_kingman_as_string="Kingman Ks sigma ({:.2E})".format(ks_predictions.theoretical_Kingman_sigma_now )
    theoretical_sigma_from_subsampling_genes_as_string = "Theoretical Ks sigma from subsampling ({:.2E})".format(
        ks_predictions.theoretical_sigma_due_to_kingman_from_subsampling_genes)
    simulated_ks_sigma_now_as_string="Simulated Ks sigma ({:.2E})".format(ks_fits.variance_from_slim)

    #plot predictions
    if show_predictions[0]:

        if config_used.DIV_time_Ge:#allopolyploid
            this_ax.plot(bins, ks_predictions.travelling_kingman_ys,
                         label='Expectation due to Kingman\nfrom Na',
                 linestyle='solid', color='b', alpha=1)
        else:
            this_ax.plot(bins, ks_predictions.autopolyploid_ys,
                         label='Expectation due to Kingman\nfrom Nb',
                 linestyle='solid', color='k',alpha=1)

    if show_predictions[1]:
        this_ax.plot(bins,ks_predictions.travelling_gaussian_ys,label='Expectation due to CLT',
                 linestyle='solid', color='r',alpha=1)

    if not include_annotation:
        return

    #add annotations
    annotation_txt = "\n".join([theoretical_ks_mean_now_as_string,
                                simulated_ks_mean_now_as_string,
                                theoretical_sigma_from_kingman_as_string,
                                theoretical_sigma_from_subsampling_genes_as_string,
                                simulated_ks_sigma_now_as_string])

    this_ax.annotate(annotation_txt, (0, 0), (0, -30), xycoords='axes fraction', textcoords='offset points', va='top')


def plot_expository_allo_auto_image(ax, png_Tnow, caption, row_idx):
    im = mpimg.imread(png_Tnow)
    ax[row_idx,0].imshow(im)
    ax[row_idx,0].get_xaxis().set_visible(False)
    ax[row_idx,0].get_yaxis().set_visible(False)
    # Selecting the axis-X making the bottom and top axes False.
    ax[row_idx,0].tick_params(axis='x', which='both', bottom=False,
                         top=False, labelbottom=False)
    ax[row_idx,0].tick_params(axis='y', which='both', right=False,
                         left=False, labelleft=False)
    for pos in ['right', 'top', 'bottom', 'left']:
        ax[row_idx,0].spines[pos].set_visible(False)
    ax[row_idx,0].set(title=caption)

def plot_expository_allo_auto_image_list(ax, images,captions):
    for i in range(0,len(images)):
        im = mpimg.imread(images[i])
        ax[i, 0].imshow(im)
        ax[i, 0].get_xaxis().set_visible(False)
        ax[i, 0].get_yaxis().set_visible(False)
        # Selecting the axis-X making the bottom and top axes False.
        ax[i, 0].tick_params(axis='x', which='both', bottom=False,
                             top=False, labelbottom=False)
        ax[i, 0].tick_params(axis='y', which='both', right=False,
                             left=False, labelleft=False)
        for pos in ['right', 'top', 'bottom', 'left']:
            ax[i, 0].spines[pos].set_visible(False)
        ax[i, 0].set(title=captions[i])

def plot_expository_allo_auto_images(ax, png_Tc, png_Tnow):


        im = mpimg.imread(png_Tnow)
        ax[0, 0].imshow(im)
        ax[0, 0].get_xaxis().set_visible(False)
        ax[0, 0].get_yaxis().set_visible(False)
        # Selecting the axis-X making the bottom and top axes False.
        ax[0, 0].tick_params(axis='x', which='both', bottom=False,
                             top=False, labelbottom=False)
        ax[0, 0].tick_params(axis='y', which='both', right=False,
                             left=False, labelleft=False)
        for pos in ['right', 'top', 'bottom', 'left']:
            ax[0, 0].spines[pos].set_visible(False)
        ax[0, 0].set(title="polyploid Ks at T_now for allo and autopolyploid\n\n\n")

        # img = Image.open(png_Tdiv)
        # im = plt.imread(get_sample_data(png_Tdiv))
        # img.close()
        im = mpimg.imread(png_Tc)
        ax[1, 0].imshow(im)
        ax[1, 0].get_xaxis().set_visible(False)
        ax[1, 0].get_yaxis().set_visible(False)
        ax[1, 0].set(title="ancestral Tc at T_div for allo and T_wgd for auto")
        for pos in ['right', 'top', 'bottom', 'left']:
            ax[1, 0].spines[pos].set_visible(False)


def add_mrca_annotations_for_autos_and_allos(this_ax, allo_config_used,
                                             avg_simulated_allo_Tc,
                                             allo_run_duration_in_m, auto_run_duration_in_m,
                                             allo_version, auto_version):
    simulated_mean_Ks_from_Tc = avg_simulated_allo_Tc * allo_config_used.Ks_per_YR
    Tc_info = 'mean Tc by Kingman = ' + \
              str(2.0 * allo_config_used.ancestral_Ne)
    ks_info = 'simulated mean Ks at Tdiv = ' + \
              "{:.2E}".format(simulated_mean_Ks_from_Tc)
    theoretical_mean_Ks_from_Tc = '2*Ne*Ks_per_YR = ' + str(allo_config_used.mean_Ks_from_Tc)
    # "{:.2E}".format(2.0 * Ne * Ks_per_YR)self.mean_Ks_from_Tc
    mut_info = 'simulated num mutations per gene = ' + \
               "{:.2E}".format(simulated_mean_Ks_from_Tc * allo_config_used.gene_length_in_bases)

    SLiM_run_time_info = ("allo run time: " + str(round(allo_run_duration_in_m, 2)) +
                          " min, version " + allo_version)

    SpecKS_run_time_info = ("auto run time: " + str(round(auto_run_duration_in_m, 2)) +
                            " min, version " + auto_version)

    annotation_txt = "\n".join([Tc_info, ks_info,
                                theoretical_mean_Ks_from_Tc, mut_info, "\n",
                                SLiM_run_time_info, SpecKS_run_time_info]) + "\n"
    this_ax.annotate(annotation_txt, (0, 0), (0, -60), xycoords='axes fraction', textcoords='offset points', va='top')


if __name__ == '__main__':
    unittest.main()

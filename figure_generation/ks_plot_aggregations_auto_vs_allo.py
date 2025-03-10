import glob
import math
import os
import unittest
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt, image as mpimg
import config
from figure_generation import ks_modeling

from figure_generation.coalescent_plot_aggregation import get_run_time_in_minutes, read_data_csv,\
    add_mrca_annotations

from figure_generation.histogram_plotter import read_Ks_csv

def plot_mrca_for_Autos_and_Allos(this_ax, allo_mrcas_by_gene, auto_mrcas_by_gene,
                                  theoretical_mrcas_by_gene,
                                  title, Ne, bin_size, xmax, ymax, total_num_genes,include_annotatons):

    num_allo_genes = len(allo_mrcas_by_gene)
    if num_allo_genes == 0:
        avg_simulated_allo_Tc = 0
    else:
        avg_simulated_allo_Tc = sum(allo_mrcas_by_gene) / num_allo_genes

    if auto_mrcas_by_gene and (len(auto_mrcas_by_gene) > 0):
        auto_mrcas_genes = len(auto_mrcas_by_gene)
        avg_simulated_auto_Tc = sum(auto_mrcas_by_gene) / auto_mrcas_genes
    else:
        avg_simulated_auto_Tc = 0

    if not xmax:
        xmax = max(allo_mrcas_by_gene)
    bins = np.arange(0, xmax, bin_size)

    two_Ne = 2.0 * Ne
    print("Tc plot bin_size_in_time: " + str(bin_size))
    kingman = [min(total_num_genes,
                   (bin_size * total_num_genes / two_Ne) * math.e ** ((-1 * i) / two_Ne))
               for i in bins]

    if allo_mrcas_by_gene:

        my_label ='Allo Tcoal by gene'
        if include_annotatons:
            my_label ='Allo Tcoal by gene\n' \
                           + "(" + str(num_allo_genes) + " genes in genome,\n" \
                           + "avg Tc " + str(int(avg_simulated_allo_Tc)) + " generations)"
        this_ax.hist(allo_mrcas_by_gene, bins=bins, facecolor='b', alpha=0.25,
                     label=my_label,
                     density=False)

    if auto_mrcas_by_gene:



        auto_mrcas_genes = len(auto_mrcas_by_gene)
        my_label ='Auto Tcoal by gene'
        if include_annotatons:
            my_label ='Auto Tcoal by gene\n' \
                           + "(" + str(auto_mrcas_genes) + " genes in genome),\n"\
                           + "avg Tc " + str(int(avg_simulated_auto_Tc)) + " generations)"
        this_ax.hist(auto_mrcas_by_gene, bins=bins, facecolor='g', alpha=0.5,
                     label=my_label ,
                     density=False)

    this_ax.plot(bins, kingman, c='red', label='Expectations under Kingman,\n'
                                               + "avg Tc " + str(int(two_Ne)) + " generations)")

    if ymax:
        this_ax.set(ylim=[0, ymax])
    this_ax.set(xlim=[0, xmax])
    this_ax.set(xlabel="MRCA time")
    this_ax.set(title=title)
    this_ax.legend()

    return avg_simulated_allo_Tc


def plot_allo_vs_auto_ks(this_ax, allo_config_used, allo_ks_by_gene,
                         auto_config_used, auto_ks_by_gene,
                         title, bin_size, xmax, ymax, show_predictions,include_annotation):

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
    theoretical_ks_mean_now= allo_config_used.mean_Ks_from_Tc + allo_config_used.t_div_as_ks
    this_ax.axvline(x=theoretical_ks_mean_now, color='purple', linestyle='--',
                    label="Expected Ks mean (Allo)",lw=3)

    this_ax.axvline(x=auto_config_used.mean_Ks_from_Nb, color='k', linestyle='--',
                    label="Expected Ks mean (Auto)",lw=3)

    #mean_Ks_from_Nb_string=  "({:.2E})".format(config_used.mean_Ks_from_Nb)
    #this_ax.axvline(x=config_used.mean_Ks_from_Nb,
    #                color='g', linestyle='--', label="Tc due to Nb " + mean_Ks_from_Nb_string )

    if include_annotation:
        add_Ks_annotations(this_ax, allo_config_used, allo_ks_by_gene, dgx_hist_ys, bins, show_predictions)

    if ymax:
        this_ax.set(ylim=[0, ymax])
    this_ax.set(xlim=[0, xmax])
    this_ax.set(xlabel="Ks")
    this_ax.set(title=title)
    this_ax.legend()


def make_Tc_Ks_Allo_vs_Auto_fig_with_subplots(num_plot_rows, bin_sizes_Ks, bin_sizes_Tc,
                                              allo_out_path, allo_run_list, run_list_name,
                                              auto_run_list, auto_out_path,
                                              xmax_Ks, xmax_Tc, ymax_Ks, ymax_Tc,
                                              suptitle, show_KS_predictions, include_annotation,
                                              plot_title_lamda):

    num_runs = len(allo_run_list)
    png_out = os.path.join(auto_out_path, run_list_name)
    #png_out = os.path.join(demographiKS_out_path, "ks_hist_by_{0}_test.jpg".format(run_list_name))
    par_dir = Path(__file__).parent.parent
    image_folder = os.path.join(par_dir, "images")
    png_Auto= os.path.join(image_folder, 'JustAutoNow.png')
    png_Tnow = os.path.join(image_folder, 'Auto_vs_Allo_now_crop_2.png')
    png_Tc = os.path.join(image_folder, 'Auto_vs_Allo_Tc_crop_2.png')

    if include_annotation:
        fig, ax = plt.subplots(num_plot_rows, num_runs, figsize=(40, 20))
        dpi_needed = 350
    else:
        fig, ax = plt.subplots(num_plot_rows, num_runs, figsize=(20,8))
        dpi_needed = 100

    fig.suptitle(suptitle)
    #captions=[
    #    "polyploid Ks at T_now for autopolyploid\n\n\n",
    #    "polyploid Ks at T_now for allo and autopolyploid\n\n\n",
    #    "ancestral Tc at T_div for allo and T_wgd for auto"]
    captions=[
        "polyploid Ks at T_now\nfor Autopolyploid",
        "polyploid Ks at T_now\nfor Allo and Autopolyploid",
        "ancestral Tc at T_div\nfor Allo and T_wgd for Auto"]
    plot_expository_allo_auto_image_list(ax, [png_Tnow,png_Tc], captions[1:3])

    joint_allo_auto_plot_index = 0
    Tc_plot_index = 1
    if num_plot_rows == 3:
        joint_allo_auto_plot_index = 1
        Tc_plot_index = 2
        plot_expository_allo_auto_image_list(ax, [png_Auto, png_Tnow, png_Tc], captions)

    for i in range(1, num_runs):
        allo_run_name = allo_run_list[i]

        if allo_run_name:

            allo_run_path = os.path.join(allo_out_path, allo_run_name)
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

            plot_title = str(plot_title_lamda(allo_config_used))
            if include_annotation:
                plot_title = "Ks at Tnow\n" + "burnin time=" + str(allo_config_used.burnin_time) + " gen,\n" \
                     + "Na=" + str(allo_config_used.ancestral_Ne)  + ", Nb=" + str(allo_config_used.bottleneck_Ne) +\
                         ", Tdiv=" + str(allo_config_used.DIV_time_Ge) + ", RC=" + str(allo_config_used.recombination_rate)


        else:
            allo_config_used = False
            allo_ks_results = []
            allo_run_duration_in_m = 0
            plot_title = "foo - didnt load a config"


        auto_run_name = auto_run_list[i]
        if auto_run_name:

            auto_run_path = os.path.join(auto_out_path,auto_run_name)
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
        else:
            spx_ks_results = []
            spx_run_duration_in_m = 0
            specks_mrcas_by_gene = False

        if num_plot_rows == 3:
            plot_allo_vs_auto_ks(ax[0, i], allo_config_used, [],
                             auto_config_used,auto_ks_results,
                             plot_title,
                             bin_sizes_Ks[0][i], xmax_Ks[0][i],
                             ymax_Ks[0][i], show_KS_predictions,
                             include_annotation)

        plot_allo_vs_auto_ks(ax[joint_allo_auto_plot_index, i], allo_config_used, allo_ks_results,
                             auto_config_used, auto_ks_results,
                plot_title, bin_sizes_Ks[joint_allo_auto_plot_index][i],
                             xmax_Ks[joint_allo_auto_plot_index][i],
                             ymax_Ks[joint_allo_auto_plot_index][i], show_KS_predictions,
                             include_annotation)


        allo_Tc_csv_file = os.path.join(allo_run_path, "1_10_simulated_ancestral_gene_mrcas.csv")
        if not os.path.exists(allo_Tc_csv_file):
            allo_Tc_csv_file = os.path.join(allo_run_path, "simulated_ancestral_gene_mrcas.csv")

        loci, allo_mrcas_by_gene = read_data_csv(allo_Tc_csv_file)

        auto_Tc_csv_file = os.path.join(auto_run_path, "1_10_simulated_ancestral_gene_mrcas.csv")
        if not os.path.exists(auto_Tc_csv_file):
            auto_Tc_csv_file = os.path.join(auto_run_path, "simulated_ancestral_gene_mrcas.csv")

        loci, auto_mrcas_by_gene = read_data_csv(auto_Tc_csv_file)

        plot_title = "Tcoal at Tdiv\nburnin time=" + str(allo_config_used.burnin_time) + " gen, " \
                     + "Na=" + str(allo_config_used.ancestral_Ne)

        theory_mrcas_by_gene=False
        avg_simulated_allo_Tc = plot_mrca_for_Autos_and_Allos(ax[Tc_plot_index, i],
                allo_mrcas_by_gene, auto_mrcas_by_gene, theory_mrcas_by_gene,
                  plot_title, allo_config_used.ancestral_Ne,
                  bin_sizes_Tc[i], xmax_Tc[i], ymax_Tc[i], allo_config_used.num_genes,
                                                              include_annotation)

        if include_annotation:
            add_mrca_annotations_for_autos_and_allos(ax[Tc_plot_index, i],
                allo_config_used, avg_simulated_allo_Tc,
                             allo_run_duration_in_m,
                             auto_run_duration_in_m,
                             allo_version, auto_version)

    ax[0, 1].set(ylabel="# paralog pairs in bin")
    ax[1, 1].set(ylabel="# genes in bin")
    plt.tight_layout()
    plt.savefig(png_out, dpi=dpi_needed)
    plt.clf()
    plt.close()

def add_Ks_annotations(this_ax, config_used,dgks_Ks_results,dgks_hist_results,
                       bins,show_predictions):

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
        this_ax.plot(bins, ks_predictions.travelling_kingman_ys, label='expectation due to Kingman',
                 linestyle='solid', color='b', alpha=1)
    if show_predictions[1]:
        this_ax.plot(bins,ks_predictions.travelling_gaussian_ys,label='expectation due to CLT',
                 linestyle='solid', color='r',alpha=1)

    #add annotations
    annotation_txt = "\n".join([theoretical_ks_mean_now_as_string,
                                simulated_ks_mean_now_as_string,
                                theoretical_sigma_from_kingman_as_string,
                                theoretical_sigma_from_subsampling_genes_as_string,
                                simulated_ks_sigma_now_as_string])

    this_ax.annotate(annotation_txt, (0, 0), (0, -30), xycoords='axes fraction', textcoords='offset points', va='top')


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

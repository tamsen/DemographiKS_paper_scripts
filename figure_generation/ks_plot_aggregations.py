import glob
import math
import os
import unittest
from pathlib import Path
import numpy as np
from PyQt5.uic.Compiler.qtproxies import i18n_func
from matplotlib import pyplot as plt
import matplotlib.image as mpimg

import config
from figure_generation import ks_modeling, curve_fitting
from figure_generation.coalescent_plot_aggregation import get_run_time_in_minutes, read_data_csv, plot_mrca, \
    add_mrca_annotations
from figure_generation.curve_fitting import wgd_lognorm2
from figure_generation.histogram_plotter import read_Ks_csv

def make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                 demographiKS_out_path, demographics_run_list, run_list_name,
                                 specks_run_list, specks_out_path,
                                 xmax_Ks, xmax_Tc, ymax_Ks, ymax_Tc,
                                 suptitle, show_KS_predictions, include_annotation, plot_title_lamda):

    num_runs = len(demographics_run_list)
    png_out = os.path.join(demographiKS_out_path, run_list_name)
    par_dir = Path(__file__).parent.parent
    image_folder = os.path.join(par_dir, "images")
    png_Tnow = os.path.join(image_folder, 'Ks_now_time_slice.jpg')
    png_Tdiv = os.path.join(image_folder, 'Tdiv_TimeSlice.jpg')
    png_with_migration_Tnow = os.path.join(image_folder, 'Migration_now.png')
    png_with_migration_Tdiv = os.path.join(image_folder, 'Migration_Tc.png')

    #fig, ax = plt.subplots(2, num_runs, figsize=(40, 20))
    num_plot_rows= 2
    if include_annotation:
        fig, ax = plt.subplots(num_plot_rows, num_runs, figsize=(40, 20))
        dpi_req = 350
    else:
        fig, ax = plt.subplots(num_plot_rows, num_runs, figsize=(20, 8))
        dpi_req = 100

    fig.suptitle(suptitle)

    for i in range(1, num_runs):
        dgx_run_name = demographics_run_list[i]

        if dgx_run_name:

            dgx_run_path = os.path.join(demographiKS_out_path, dgx_run_name)
            print("dgx_run_path: " +dgx_run_path )
            glob_results=glob.glob(dgx_run_path + '/*.used.xml')
            input_xml_file = glob_results[0]
            config_used = config.DemographiKS_config(input_xml_file)
            csv_file_name = config_used.sim_name + '.csv'
            ks_file = os.path.join(dgx_run_path, csv_file_name)
            print("reading " + ks_file)
            demographiKS_ks_results = read_Ks_csv(ks_file, False)
            dgx_run_duration_in_m, dgx_version = get_run_time_in_minutes(dgx_run_path)

            if include_annotation:
                plot_title = "Ks at Tnow\n" + \
                         "burnin time=" + str(config_used.burnin_time) + " gen, " \
                     + "Na=" + str(config_used.ancestral_Ne)  + ", Nb=" + str(config_used.bottleneck_Ne) +\
                         ",\nTdiv=" + str(config_used.DIV_time_Ge) + ", RC=" + str(config_used.recombination_rate)
                if config_used.mig_rate:
                    plot_title = plot_title + ", MigStart=" + str(config_used.mig_start) + \
                              ", MigEnd=" + str(config_used.mig_stop) +  ", MigRate=" + str(config_used.mig_rate)
                else:
                    plot_title = plot_title + ", MigRate=0"
            else:
                plot_title = str(plot_title_lamda(config_used))

        else:
            config_used = False
            demographiKS_ks_results = []
            dgx_run_duration_in_m = 0
            plot_title = "Ks at Tnow"
            dgx_version = "NA"

        spx_run_name = specks_run_list[i]
        if spx_run_name:
            spx_run_nickname = spx_run_name.split('_')[1]
            spx_run_path = os.path.join(specks_out_path, spx_run_name)
            csv_file_name = 'Allo_' + spx_run_nickname + '_ML_rep0_Ks_by_GeneTree.csv'
            spx_ks_results = read_Ks_csv(os.path.join(spx_run_path,csv_file_name), True)
            spx_run_duration_in_m,spx_version = get_run_time_in_minutes(spx_run_path)
            specks_csv_file = os.path.join(spx_run_path, "variations_in_div_time.txt")
            loci, specks_mrcas_by_gene = read_data_csv(specks_csv_file)

            glob_results=glob.glob(spx_run_path + '/*.used.xml')
            input_xml_file = glob_results[0]
            specks_config_used = config.DemographiKS_config(input_xml_file)

        else:
            spx_ks_results = []
            spx_run_duration_in_m = 0
            spx_version= "NA"
            specks_mrcas_by_gene = False


        plot_ks(ax[0, i], config_used, demographiKS_ks_results, spx_ks_results,
                plot_title, bin_sizes_Ks[i], xmax_Ks[i], ymax_Ks[i],
                show_KS_predictions, include_annotation)


        if dgx_run_name:
            slim_csv_file_0 = os.path.join(dgx_run_path, "simulated_ancestral_gene_mrcas.csv")
            if os.path.exists(slim_csv_file_0):
                loci, slim_mrcas_by_gene = read_data_csv(slim_csv_file_0)
            else:
                slim_csv_file_1 = os.path.join(dgx_run_path, "1_5_simulated_ancestral_gene_mrcas.csv")
                loci, slim_mrcas_by_gene = read_data_csv(slim_csv_file_1)
        else:
            slim_mrcas_by_gene=[]

        #theory_output_file = os.path.join(dgx_run_path, "theoretical_ancestral_gene_mrcas.csv")
        #loci, theory_mrcas_by_gene = read_data_csv(theory_output_file)
        if include_annotation:
            plot_title = "Tcoal at Tdiv\nburnin time=" + str(config_used.burnin_time) + " gen, " \
                     + "Na=" + str(config_used.ancestral_Ne)
        else:
            plot_title = "Tcoal at Tdiv"

        theory_mrcas_by_gene=False
        avg_slim_Tc = plot_mrca(ax[1, i], slim_mrcas_by_gene, specks_mrcas_by_gene, theory_mrcas_by_gene,
                  plot_title, config_used.ancestral_Ne,
                  bin_sizes_Tc[i], xmax_Tc[i], ymax_Tc[i], config_used.num_genes, include_annotation)

        if include_annotation:
            add_mrca_annotations(ax[1, i], config_used, avg_slim_Tc,
                             dgx_run_duration_in_m, spx_run_duration_in_m,
                             dgx_version,spx_version)

    if config_used.mig_rate > 0:
        plot_expository_images(ax, png_with_migration_Tdiv , png_with_migration_Tnow)
    else:
        plot_expository_images(ax, png_Tdiv, png_Tnow)

    ax[0, 1].set(ylabel="# paralog pairs in bin")
    ax[1, 1].set(ylabel="# genes in bin")
    plt.tight_layout()
    plt.savefig(png_out, dpi=dpi_req)
    plt.clf()
    plt.close()


def plot_expository_images(ax, png_Tdiv, png_Tnow):

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
        ax[0, 0].set(title="polyploid Ks at T_now")

        #img = Image.open(png_Tdiv)
        #im = plt.imread(get_sample_data(png_Tdiv))
        #img.close()
        im = mpimg.imread(png_Tdiv)
        ax[1, 0].imshow(im)
        ax[1, 0].get_xaxis().set_visible(False)
        ax[1, 0].get_yaxis().set_visible(False)
        ax[1, 0].set(title="ancestral Tc at T_div")
        for pos in ['right', 'top', 'bottom', 'left']:
            ax[1, 0].spines[pos].set_visible(False)


def plot_ks(this_ax, config_used, slim_ks_by_gene, spx_ks_by_gene,
            title, bin_size, xmax, ymax, show_predictions, include_annotations):

    num_slim_genes = len(slim_ks_by_gene)
    num_specks_genes = len(spx_ks_by_gene)
    include_logfit=False
    include_RC_model=False

    if not xmax:
        xmax = max(slim_ks_by_gene)
    bins = np.arange(0, xmax, bin_size)

    if len(slim_ks_by_gene) > 0:
        if include_annotations:
            my_label='SLiM Ks by gene'
        else:
            my_label='SLiM Ks by gene\n'\
                           + "(" + str(num_slim_genes) + " paralogs in genome)"
        dgx_hist_ys, bins, patches = this_ax.hist(slim_ks_by_gene, bins=bins, facecolor='b', alpha=0.25,
                     label=my_label,
                     density=False)

    if len(spx_ks_by_gene) > 0:

        if include_annotations:
            my_label='SpecKS Ks by gene'
        else:
            my_label='SpecKS Ks by gene\n'\
                           + "(" + str(num_specks_genes) + " paralogs in genome)"
        this_ax.hist(spx_ks_by_gene, bins=bins, facecolor='c', alpha=0.25,
                     label=my_label,
                     density=False)

    half_bin_size=0.5*bin_size
    if config_used.t_div_as_ks:
        this_ax.axvline(x=half_bin_size+config_used.t_div_as_ks, color='b', linestyle='-.', label="input Tdiv as Ks")

    if not config_used.DIV_time_Ge: #this would be an autopolyploid
        this_ax.axvline(x=half_bin_size + config_used.mean_Ks_from_Nb, color='g', linestyle='--',
                        label="Expected Ks mean (4Nb)")
    theoretical_ks_mean_now=config_used.mean_Ks_from_Tc+config_used.t_div_as_ks
    this_ax.axvline(x=half_bin_size+theoretical_ks_mean_now, color='r', linestyle='--',
                    label="Expected Ks mean (2Na)")

    RC_rate_per_base_per_year=config_used.recombination_rate
    exp_num_RC_per_base_since_DIV=RC_rate_per_base_per_year*theoretical_ks_mean_now/config_used.Ks_per_YR
    exp_num_RC_per_gene = exp_num_RC_per_base_since_DIV * config_used.gene_length_in_bases
    exp_num_RC_per_gene_as_int = round(exp_num_RC_per_gene,4)
    max_x_exp_using_RC = config_used.t_div_as_ks + (config_used.mean_Ks_from_Tc * exp_num_RC_per_gene / (exp_num_RC_per_gene + 1))

    if include_RC_model:
            this_ax.axvline(x=half_bin_size+max_x_exp_using_RC, color='g', linestyle=':',
                    label="Predicted Ks peak for\nRC=" + str(exp_num_RC_per_gene_as_int))
            print("predicted KS peak:\t" + str(half_bin_size+max_x_exp_using_RC))

    if len(slim_ks_by_gene) > 0:
        dgx_hist_ys_list=list(dgx_hist_ys)
        max_value = max(dgx_hist_ys_list)
        max_index = dgx_hist_ys_list.index(max_value)
        max_x = bins[max_index]
        print("true KS peak:\t" + str(max_x))
        print("half_bin_size:\t" + str(half_bin_size))

    # a good start..
    if len(slim_ks_by_gene) > 0:
        shape = 0.5
        num_genes= num_slim_genes
        loc = -1*config_used.ancestral_Ne*config_used.Ks_per_YR #-1000  # should be -N
        scale = config_used.ancestral_Ne * math.sqrt(2.0*math.pi)*config_used.Ks_per_YR
        p0=[bin_size*num_genes,shape,loc,scale]
        fit_curve_ys_ln2, xs_for_wgd, popt = \
            curve_fitting.fit_curve_to_xs_and_ys(bins[0:-1],
                                             dgx_hist_ys, curve_fitting.wgd_lognorm2, p0=p0)
        if include_logfit:
            if fit_curve_ys_ln2:
                this_ax.plot(xs_for_wgd,fit_curve_ys_ln2,color='g', label="Ks logfit")

    else:
        fit_curve_ys_ln2 = False
        dgx_hist_ys = False



    #mean_Ks_from_Nb_string=  "({:.2E})".format(config_used.mean_Ks_from_Nb)
    #this_ax.axvline(x=config_used.mean_Ks_from_Nb,
    #                color='g', linestyle='--', label="Tc due to Nb " + mean_Ks_from_Nb_string )
    if config_used.mig_rate > 0:
        mig_start_as_Ks=config_used.t_div_as_ks- (config_used.Ks_per_YR*config_used.mig_start)
        mig_stop_as_Ks=config_used.t_div_as_ks- (config_used.Ks_per_YR*config_used.mig_stop)
        this_ax.axvline(x=mig_start_as_Ks, color='g', linestyle=':', label="Mig start")
        this_ax.axvline(x=mig_stop_as_Ks, color='g', linestyle='--', label="Mig stop")
        this_ax.axvspan(mig_stop_as_Ks, mig_start_as_Ks, alpha=0.25, color='g')

    if len(slim_ks_by_gene) > 0:
        if include_annotations:
            add_Ks_annotations(this_ax, config_used, slim_ks_by_gene,dgx_hist_ys, bins, show_predictions)

    if ymax:
        this_ax.set(ylim=[0, ymax])
    this_ax.set(xlim=[0, xmax])
    this_ax.set(xlabel="Ks")
    this_ax.set(title=title)
    this_ax.legend()


def add_Ks_annotations(this_ax, config_used,dgks_Ks_results,dgks_hist_results,
                       bins,show_predictions):

    ks_predictions = ks_modeling.Ks_modeling_predictions(config_used, bins)
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


if __name__ == '__main__':
    unittest.main()

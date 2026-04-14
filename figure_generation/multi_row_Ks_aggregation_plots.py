import glob
import math
import os
import unittest
import numpy as np
from matplotlib import pyplot as plt

import config
from figure_generation import ks_modeling, curve_fitting
from figure_generation.coalescent_plot_aggregation import get_run_time_in_minutes, read_data_csv
from figure_generation.histogram_plotter import read_Ks_csv

class SPKSvsDGKSPlotData:
        def __init__(self, demographiKS_out_path, demographics_run_list,demographics_run_list_name,
                     specks_out_path,specks_run_list,
                     bin_sizes_Ks,
                     xmax_Ks, ymax_Ks,
                     show_KS_predictions,
                     include_annotation, which_plot_panels_to_show_legend, plot_title_lamda):
            self.demographiKS_out_path = demographiKS_out_path
            self.demographics_run_list = demographics_run_list
            self.demographics_run_list_name = demographics_run_list_name
            self.specks_out_path = specks_out_path
            self.specks_run_list = specks_run_list
            self.bin_sizes_Ks = bin_sizes_Ks
            self.xmax_Ks=xmax_Ks
            self.ymax_Ks=ymax_Ks
            self.show_KS_predictions=show_KS_predictions
            self.include_annotation=include_annotation
            self.which_plot_panels_to_show_legend=which_plot_panels_to_show_legend
            self.plot_title_lamda=plot_title_lamda

def make_multi_row_Ks_fig_with_subplots(plotdatalist, output_png_path):

    first_plot = plotdatalist[0]
    num_rows=len(plotdatalist)
    num_runs = len(first_plot.auto_run_list)

    fig, ax = plt.subplots(num_rows, num_runs, figsize=(20, 16))
    dpi_req = 100

    for r in range(0, num_rows):
        this_plot_data = plotdatalist[r]
        for i in range(0, num_runs):
            dgx_run_name = this_plot_data.auto_run_list[i]

            if dgx_run_name:

                dgx_run_path = os.path.join(this_plot_data.auto_out_path, dgx_run_name)
                print("dgx_run_path: " +dgx_run_path )
                glob_results=glob.glob(dgx_run_path + '/*.used.xml')
                input_xml_file = glob_results[0]
                config_used = config.DemographiKS_config(input_xml_file)
                csv_file_name = config_used.sim_name + '.csv'
                ks_file = os.path.join(dgx_run_path, csv_file_name)
                print("reading " + ks_file)
                demographiKS_ks_results = read_Ks_csv(ks_file, False)
                dgx_run_duration_in_m, dgx_version = get_run_time_in_minutes(dgx_run_path)

                if this_plot_data.include_annotation:
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
                    plot_title = str(this_plot_data.plot_title_lamda(config_used))


            else:
                config_used = False
                demographiKS_ks_results = []
                dgx_run_duration_in_m = 0
                plot_title = "Ks at Tnow"
                dgx_version = "NA"

            spx_run_name = this_plot_data.allo_run_list[i]
            if spx_run_name:
                spx_run_nickname = spx_run_name.split('_')[1]
                spx_run_path = os.path.join(this_plot_data.allo_out_path, spx_run_name)
                csv_file_name = 'Allo_' + spx_run_nickname + '_ML_rep0_Ks_by_GeneTree.csv'
                spx_ks_results = read_Ks_csv(os.path.join(spx_run_path,csv_file_name), True)
                spx_run_duration_in_m,spx_version = get_run_time_in_minutes(spx_run_path)
                specks_csv_file = os.path.join(spx_run_path, "variations_in_div_time.txt")
                loci, specks_mrcas_by_gene = read_data_csv(specks_csv_file)
                print(spx_run_path)
                glob_results=glob.glob(spx_run_path + '/*.used.xml')
                input_xml_file = glob_results[0]
                #specks_config_used = config.DemographiKS_config(input_xml_file)
                if not config_used:
                    config_used = config.DemographiKS_config(input_xml_file)

                if not dgx_run_name:
                    if this_plot_data.include_annotation:
                        plot_title = "Ks at Tnow\n"
                    else:
                        plot_title = str(this_plot_data.plot_title_lamda(config_used))

            else:
                spx_ks_results = []



            dgx_hist_ys, bins = plot_ks(i,ax[r, i], config_used, demographiKS_ks_results, spx_ks_results,
                    plot_title, this_plot_data.bin_sizes_Ks[i], this_plot_data.xmax_Ks[i],
                    this_plot_data.ymax_Ks[i],
                    this_plot_data.show_KS_predictions, this_plot_data.include_annotation,
                    this_plot_data.which_plot_panels_to_show_legend)

            #if dgx_hist_ys:
            csv_out = os.path.join(this_plot_data.auto_out_path,
                                   this_plot_data.auto_run_list_name + "_out.csv")
            with open(csv_out, 'a') as f:

                    run_name = this_plot_data.plot_title_lamda(config_used).replace("Ks at Tnow\n", "")
                    dgx_hist_ys_string = " ".join([str(d) for d in dgx_hist_ys])
                    bins_string = " ".join([str(b) for b in bins])
                    data = [run_name, bins_string, dgx_hist_ys_string]
                    f.writelines(",".join(data) + "\n")

    ax[r,0].set(ylabel="# paralog pairs in bin")

    plt.tight_layout()
    if this_plot_data.include_annotation:
        plt.savefig(output_png_path +"_annotated.png", dpi=dpi_req)
    else:
        plt.savefig(output_png_path, dpi=dpi_req)
    plt.clf()
    plt.close()



def plot_ks(i, this_ax, config_used, slim_ks_by_gene, spx_ks_by_gene,
            title, bin_size, xmax, ymax, show_predictions, include_annotations,plots_to_show_legend):

    num_slim_genes = len(slim_ks_by_gene)
    num_specks_genes = len(spx_ks_by_gene)
    include_logfit=("RC" in title)
    include_RC_model=("RC" in title)

    #if not xmax:
    #    xmax = max(slim_ks_by_gene)
    bins = np.arange(0, xmax, bin_size)

    if len(slim_ks_by_gene) > 0:
        if include_annotations:
            my_label='DGKS Ks by gene' + "(" + str(num_slim_genes) + " paralogs in genome)"
        else:
            my_label='DGKS Ks by gene'
        dgx_hist_ys, bins, patches = this_ax.hist(slim_ks_by_gene, bins=bins, facecolor='b', alpha=0.25,
                     label=my_label,
                     density=False)

    if len(spx_ks_by_gene) > 0:

        if include_annotations:
            my_label='SpecKS Ks by gene'\
                + "(" + str(num_specks_genes) + " paralogs in genome)"
        else:
            my_label='SpecKS Ks by gene'
        dgx_hist_ys, bins, patches =this_ax.hist(spx_ks_by_gene, bins=bins, facecolor='c', alpha=0.25,
                     label=my_label,
                     density=False)

    half_bin_size=0.5*bin_size
    if config_used.t_div_as_ks:
        this_ax.axvline(x=half_bin_size+config_used.t_div_as_ks, color='b', linestyle='-.',
                        label="input Tdiv as Ks")

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
            rc_label="Predicted Ks peak for\nRC=" + str(exp_num_RC_per_gene_as_int)
            this_ax.axvline(x=half_bin_size+max_x_exp_using_RC, color='g', linestyle=':',
                    label="Predicted Ks peak",lw=4)
            print("predicted KS peak:\t" + str(half_bin_size+max_x_exp_using_RC))

    if len(slim_ks_by_gene) > 0:
        dgx_hist_ys_list=list(dgx_hist_ys)
        max_value = max(dgx_hist_ys_list)
        max_index = dgx_hist_ys_list.index(max_value)
        max_x = bins[max_index]
        print("true KS peak:\t" + str(max_x))
        print("half_bin_size:\t" + str(half_bin_size))

    # a good start..
    #if len(slim_ks_by_gene) > 0:
    if include_logfit:
        bin_centers=[0.5*(bins[k] + bins[k+1]) for k in range(0,len(bins)-1)]
        shape = 0.5
        num_genes= num_slim_genes
        loc = -1*config_used.ancestral_Ne*config_used.Ks_per_YR #-1000  # should be -N
        scale = config_used.ancestral_Ne * math.sqrt(2.0*math.pi)*config_used.Ks_per_YR
        p0=[bin_size*num_genes,shape,loc,scale]
        fit_curve_ys_ln2, xs_for_wgd, popt = \
            curve_fitting.fit_curve_to_xs_and_ys(bins[0:-1],
                                             dgx_hist_ys, curve_fitting.wgd_lognorm2, p0=p0)
        if fit_curve_ys_ln2:
                this_ax.plot(bin_centers,fit_curve_ys_ln2,color='g', label="Ks logfit")
        else:
            subset_xs=[]
            subset_ys=[]

            #as a rescue, just fit around the main peak
            for j in range(0,len(dgx_hist_ys)):
                x=bins[j]
                if x>config_used.theoretical_ks_mean_now - 2*bin_size:
                    subset_xs.append(x)
                    subset_ys.append(dgx_hist_ys[j])

            fit_curve_ys_subset, xs_subset, popt = \
                curve_fitting.fit_curve_to_xs_and_ys(subset_xs,
                                                     subset_ys, curve_fitting.wgd_lognorm2, p0=p0)
            if fit_curve_ys_subset:
                fit_curve_ys_ln2 = [curve_fitting.wgd_lognorm2(x, *popt) for x in xs_for_wgd]
                this_ax.plot(bin_centers, fit_curve_ys_ln2, color='g', label="Ks lognorm fit")
                fit_curve_ys_ln2 = True
                dgx_hist_ys = True
    #else:
    #    #fit_curve_ys_ln2 = False
    #    #dgx_hist_ys = False



    #mean_Ks_from_Nb_string=  "({:.2E})".format(config_used.mean_Ks_from_Nb)
    #this_ax.axvline(x=config_used.mean_Ks_from_Nb,
    #                color='g', linestyle='--', label="Tc due to Nb " + mean_Ks_from_Nb_string )
    if config_used.mig_rate > 0:
        mig_start_as_Ks=config_used.t_div_as_ks- (config_used.Ks_per_YR*config_used.mig_start)
        mig_stop_as_Ks=config_used.t_div_as_ks- (config_used.Ks_per_YR*config_used.mig_stop)
        this_ax.axvline(x=mig_start_as_Ks, color='g', linestyle=':', label="Mig start", linewidth=3)
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
    #this_ax.legend(loc='upper left', bbox_to_anchor=(-1, 2))


    if i in plots_to_show_legend:
        this_ax.legend()

    #    if include_RC_model:
    #        this_ax.legend(loc='upper center', bbox_to_anchor=(0.5, 2.0), ncol=1)
    #    else:
    #        this_ax.legend()

    return dgx_hist_ys, bins

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

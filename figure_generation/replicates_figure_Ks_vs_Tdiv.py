import glob
import math
import os.path
import unittest
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import config
from figure_generation.coalescent_plot_aggregation \
    import get_run_time_in_minutes, read_data_csv, add_mrca_annotations

from figure_generation.histogram_plotter import read_Ks_csv
from figure_generation.ks_plot_aggregations import plot_ks, plot_expository_images

class TestResampleTc(unittest.TestCase):

    #        <DIV_time_Ge>10000</DIV_time_Ge>
    #< DIV_time_Ge > 100000 < / DIV_time_Ge >
    #< DIV_time_Ge > 500000 < / DIV_time_Ge >
    #< DIV_time_Ge > 1000000 < / DIV_time_Ge >
    def test_Replicates_With_Tdiv(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPCKS_vs_DGKS_replicates/Ne10'
        specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPCKS_vs_DGKS_replicates/Ne10'
        #/ usr / scratch2 / userdata2 / tdunn / DemographiKS_output / TE
        run_list = [
        'TE05fix_rep1_m04d12y2025_h11m07s35',
        'TE05fix_rep2_m04d12y2025_h11m16s05',
        'TE05fix_rep3_m04d12y2025_h11m16s06',
        'TE05fix_rep4_m04d12y2025_h11m16s08',
        'TE05fix_rep5_m04d12y2025_h11m16s10',
        'TE05fix_rep6_m04d12y2025_h11m16s13',
        'TE05fix_rep7_m04d12y2025_h11m16s15',
        'TE05fix_rep8_m04d12y2025_h11m16s17',
        'TE05fix_rep9_m04d12y2025_h11m16s20',
        'TE05fix_rep10_m04d12y2025_h11m16s23',
        'TE07_fix_rep1_m04d14y2025_h10m26s12',
        'TE07_fix_rep2_m04d14y2025_h10m35s02',
        'TE07_fix_rep3_m04d14y2025_h10m35s04',
        'TE07_fix_rep4_m04d14y2025_h10m35s05',
        'TE07_fix_rep5_m04d14y2025_h10m35s09',
        'TE07_fix_rep6_m04d14y2025_h10m35s10',
        'TE07_fix_rep7_m04d14y2025_h10m35s12',
        'TE07_fix_rep8_m04d14y2025_h10m35s14',
        'TE07_fix_rep9_m04d14y2025_h10m35s16',
        'TE07_fix_rep10_m04d14y2025_h10m35s19',
        'TE08_fix_rep1_m04d15y2025_h09m17s04',
        'TE08_fix_rep1_m04d15y2025_h09m21s40',
        'TE08_fix_rep1_m04d15y2025_h11m17s19',
        'TE08_fix_rep2_m04d15y2025_h11m09s33',
        'TE08_fix_rep3_m04d15y2025_h11m10s14',
        'TE08_fix_rep4_m04d15y2025_h11m16s54',
        'TE08_fix_rep5_m04d15y2025_h11m16s56',
        'TE08_fix_rep6_m04d15y2025_h11m16s58',
        'TE08_fix_rep7_m04d15y2025_h11m17s00',
        'TE08_fix_rep8_m04d15y2025_h11m17s03',
        'TE08_fix_rep9_m04d15y2025_h11m17s08',
        'TE08_fix_rep10_m04d15y2025_h11m17s12',
        'TE09_fix_rep1_m04d16y2025_h10m20s41',
        'TE09_fix_rep2_m04d16y2025_h10m20s44',
        'TE09_fix_rep3_m04d16y2025_h10m20s46',
        'TE09_fix_rep4_m04d16y2025_h10m20s48',
        'TE09_fix_rep5_m04d16y2025_h10m20s50',
        'TE09_fix_rep6_m04d16y2025_h10m20s54',
        'TE09_fix_rep7_m04d16y2025_h10m20s56',
        'TE09_fix_rep8_m04d16y2025_h10m20s58',
        'TE09_fix_rep9_m04d16y2025_h10m21s05',
        'TE09_fix_rep10_m04d16y2025_h10m21s08']

        xmax_Ks = [0.15 for f in run_list]
        bin_sizes_Ks = [xmax_Ks_i / 25 for xmax_Ks_i in xmax_Ks]
        xmax_Tc = [1000 for f in run_list ]
        bin_sizes_Tc =[xmax_Tc_i / 25 for xmax_Tc_i in xmax_Tc]
        ymax_Tc = [False for f in run_list]
        run_list_num = "DGKS_1000_gen_by_Ne_fast_mut_rate_Replicates_Ne10."
        ymax_Ks = [400 for f in run_list ]
        specks_TE_run_list = [False for f in run_list ]
        suptitle = "SLiM vs SpecKS, Tcoal and Ks"
        show_KS_predictions=[False,False,False]
        include_annotation=False
        plot_title_lamda = lambda config: "Ks at Tnow\n"+ "Ne:" + str(config.ancestral_Ne)
        which_plot_panels_to_show_legend = [1]

        make_Tc_Ks_fig_for_replicates(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, run_list, run_list_num,
                                     specks_TE_run_list, specks_out_path,
                                     xmax_Ks, xmax_Tc, ymax_Ks, ymax_Tc,
                                      suptitle, show_KS_predictions,
                                     include_annotation,which_plot_panels_to_show_legend,plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here
        plt.close()


    def test_Replicates_With_Ne100(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPCKS_vs_DGKS_replicates/Ne100'
        specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPCKS_vs_DGKS_replicates/Ne100'

        run_list = [
            'DGKS_Ne100_rep1_m03d26y2025_h15m02s34',
            'DGKS_Ne100_rep2_m03d26y2025_h15m02s36',
            'DGKS_Ne100_rep3_m03d26y2025_h15m02s38',
            'DGKS_Ne100_rep4_m03d26y2025_h15m02s39',
            'DGKS_Ne100_rep5_m03d26y2025_h15m02s41',
            'DGKS_Ne100_rep6_m03d26y2025_h15m02s43',
            'DGKS_Ne100_rep7_m03d26y2025_h15m02s45',
            'DGKS_Ne100_rep8_m03d26y2025_h15m02s47',
            'DGKS_Ne100_rep9_m03d26y2025_h15m02s49',
            'DGKS_Ne100_rep10_m03d26y2025_h15m02s52']

        xmax_Ks = [0.15 for f in run_list]
        bin_sizes_Ks = [xmax_Ks_i / 25 for xmax_Ks_i in xmax_Ks]
        xmax_Tc = [2000 for f in run_list ]
        bin_sizes_Tc =[xmax_Tc_i / 25 for xmax_Tc_i in xmax_Tc]
        ymax_Tc = [False for f in run_list]
        run_list_num = "DGKS_1000_gen_by_Ne_fast_mut_rate_Replicates_Ne100."
        ymax_Ks = [400 for f in run_list ]
        specks_TE_run_list = [False for f in run_list ]
        suptitle = "SLiM vs SpecKS, Tcoal and Ks"
        show_KS_predictions=[False,False,False]
        include_annotation=False
        plot_title_lamda = lambda config: "Ks at Tnow\n"+ "Ne:" + str(config.ancestral_Ne)
        which_plot_panels_to_show_legend = [1]

        make_Tc_Ks_fig_for_replicates(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, run_list, run_list_num,
                                     specks_TE_run_list, specks_out_path,
                                     xmax_Ks, xmax_Tc, ymax_Ks, ymax_Tc,
                                      suptitle, show_KS_predictions,
                                     include_annotation,which_plot_panels_to_show_legend,plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here
        plt.close()

    def test_Replicates_With_Ne500(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPCKS_vs_DGKS_replicates/Ne500'
        specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPCKS_vs_DGKS_replicates/Ne500'

        run_list = [
            'DGKS_Ne500_rep1_m03d26y2025_h17m09s08',
            'DGKS_Ne500_rep2_m03d26y2025_h17m09s07',
            'DGKS_Ne500_rep3_m03d26y2025_h17m09s11',
            'DGKS_Ne500_rep4_m03d26y2025_h17m09s16',
            'DGKS_Ne500_rep5_m03d26y2025_h17m09s19',
            'DGKS_Ne500_rep6_m03d26y2025_h17m09s21',
            'DGKS_Ne500_rep7_m03d26y2025_h17m09s23',
            'DGKS_Ne500_rep8_m03d26y2025_h17m09s25',
            'DGKS_Ne500_rep9_m03d26y2025_h17m09s29',
            'DGKS_Ne500_rep10_m03d26y2025_h17m09s32']

        xmax_Ks = [0.15 for f in run_list]
        bin_sizes_Ks = [xmax_Ks_i / 25 for xmax_Ks_i in xmax_Ks]
        xmax_Tc = [10000 for f in run_list ]
        bin_sizes_Tc =[xmax_Tc_i / 25 for xmax_Tc_i in xmax_Tc]
        ymax_Tc = [False for f in run_list]
        run_list_num = "DGKS_1000_gen_by_Ne_fast_mut_rate_Replicates_Ne500."
        ymax_Ks = [400 for f in run_list ]
        specks_TE_run_list = [False for f in run_list ]
        suptitle = "SLiM vs SpecKS, Tcoal and Ks"
        show_KS_predictions=[False,False,False]
        include_annotation=False
        plot_title_lamda = lambda config: "Ks at Tnow\n"+ "Ne:" + str(config.ancestral_Ne)
        which_plot_panels_to_show_legend = [1]

        make_Tc_Ks_fig_for_replicates(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, run_list, run_list_num,
                                     specks_TE_run_list, specks_out_path,
                                     xmax_Ks, xmax_Tc, ymax_Ks, ymax_Tc,
                                      suptitle, show_KS_predictions,
                                     include_annotation,which_plot_panels_to_show_legend,plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here
        plt.close()


    def test_Replicates_With_Ne1000(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPCKS_vs_DGKS_replicates/Ne1000'
        specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPCKS_vs_DGKS_replicates/Ne1000'

        run_list = [
                'DGKS_Ne1000_rep1_m04d08y2025_h13m41s15',
                'DGKS_Ne1000_rep2_m04d08y2025_h13m41s16',
                'DGKS_Ne1000_rep3_m04d08y2025_h13m41s18',
                'DGKS_Ne1000_rep4_m04d08y2025_h13m41s20',
                'DGKS_Ne1000_rep5_m04d08y2025_h13m41s22',
                'DGKS_Ne1000_rep6_m04d08y2025_h13m41s25',
                'DGKS_Ne1000_rep7_m04d08y2025_h13m41s27',
                'DGKS_Ne1000_rep8_m04d08y2025_h13m41s30',
                'DGKS_Ne1000_rep9_m04d08y2025_h13m41s32',
                'DGKS_Ne1000_rep10_m04d08y2025_h13m41s12']

        xmax_Ks = [0.15 for f in run_list]
        bin_sizes_Ks = [xmax_Ks_i / 25 for xmax_Ks_i in xmax_Ks]
        xmax_Tc = [80000 for f in run_list ]
        bin_sizes_Tc =[xmax_Tc_i / 25 for xmax_Tc_i in xmax_Tc]
        ymax_Tc = [False for f in run_list]
        run_list_num = "DGKS_1000_gen_by_Ne_fast_mut_rate_Replicates_Ne1000."
        ymax_Ks = [400 for f in run_list ]
        specks_TE_run_list = [False for f in run_list ]
        suptitle = "SLiM vs SpecKS, Tcoal and Ks"
        show_KS_predictions=[False,False,False]
        include_annotation=False
        plot_title_lamda = lambda config: "Ks at Tnow\n"+ "Ne:" + str(config.ancestral_Ne)
        which_plot_panels_to_show_legend = [1]

        make_Tc_Ks_fig_for_replicates(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, run_list, run_list_num,
                                     specks_TE_run_list, specks_out_path,
                                     xmax_Ks, xmax_Tc, ymax_Ks, ymax_Tc,
                                      suptitle, show_KS_predictions,
                                     include_annotation,which_plot_panels_to_show_legend,plot_title_lamda)

        self.assertEqual(True, True)
        plt.close()

def make_Tc_Ks_fig_for_replicates(bin_sizes_Ks, bin_sizes_Tc,
                                 demographiKS_out_path, demographics_run_list, run_list_name,
                                 specks_run_list, specks_out_path,
                                 xmax_Ks, xmax_Tc, ymax_Ks, ymax_Tc,
                                 suptitle, show_KS_predictions, include_annotation,
                                 plots_to_show_legend,
                                 plot_title_lamda):

    num_runs = len(demographics_run_list)
    png_out = os.path.join(demographiKS_out_path, run_list_name)
    csv_out = png_out+"csv"
    par_dir = Path(__file__).parent.parent
    image_folder = os.path.join(par_dir, "images")
    png_Tnow = os.path.join(image_folder, 'Ks_now_time_slice.jpg')
    png_Tdiv = os.path.join(image_folder, 'Tdiv_TimeSlice.jpg')
    png_with_migration_Tnow = os.path.join(image_folder, 'Ks_with_migration.png')
    png_with_migration_Tdiv = os.path.join(image_folder, 'Tc_with_migration.png')

    num_rows=2
    if include_annotation:
        fig, ax = plt.subplots(num_rows, 2, figsize=(20, 20))
        dpi_req = 350
    else:
        fig, ax = plt.subplots(num_rows, 2, figsize=(10, 8))
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
            #specks_config_used = config.DemographiKS_config(input_xml_file)
            if not config_used:
                config_used = config.DemographiKS_config(input_xml_file)

            if not dgx_run_name:
                if include_annotation:
                    plot_title = "Ks at Tnow\n"
                else:
                    plot_title = str(plot_title_lamda(config_used))

        else:
            spx_ks_results = []
            spx_run_duration_in_m = 0
            spx_version= "NA"
            specks_mrcas_by_gene = False


        dgx_hist_ys, bins=plot_ks(i,ax[0, 1], config_used, demographiKS_ks_results, spx_ks_results,
                plot_title, bin_sizes_Ks[i], xmax_Ks[i], ymax_Ks[i],
                show_KS_predictions, include_annotation,plots_to_show_legend)

        with open(csv_out, 'a') as f:

            run_name=plot_title_lamda(config_used).replace("Ks at Tnow\n","")
            if include_annotation:
                dgx_hist_ys_string=" ".join( [str(d) for d in dgx_hist_ys])
                bins_string=" ".join([str(b) for b in bins])
                data = [run_name, bins_string, dgx_hist_ys_string]
                f.writelines(",".join(data) + "\n")

        if num_rows < 2:
            continue

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
        avg_slim_Tc = plot_replicate_mrca(ax[1, 1], slim_mrcas_by_gene, specks_mrcas_by_gene, theory_mrcas_by_gene,
                  plot_title, config_used.ancestral_Ne,
                  bin_sizes_Tc[i], xmax_Tc[i], ymax_Tc[i], config_used.num_genes, include_annotation)

        if include_annotation:
            add_mrca_annotations(ax[1, i], config_used, avg_slim_Tc,
                             dgx_run_duration_in_m, spx_run_duration_in_m,
                             dgx_version,spx_version)

    if config_used.mig_rate > 0:
        plot_expository_images(num_rows, ax, png_with_migration_Tdiv , png_with_migration_Tnow)
    else:
        plot_expository_images(num_rows, ax, png_Tdiv, png_Tnow)

    ax[0, 1].set(ylabel="# paralog pairs in bin")

    if num_rows > 1:
        ax[1, 1].set(ylabel="# genes in bin")
    plt.tight_layout()
    if include_annotation:
        plt.savefig(png_out +"_annotated.png", dpi=dpi_req)
    else:
        plt.savefig(png_out, dpi=dpi_req)
    plt.clf()
    plt.close()


def plot_replicate_mrca(this_ax, slim_mrcas_by_gene, specks_mrcas_by_gene, theoretical_mrcas_by_gene,
              title, Ne, bin_size, xmax, ymax, total_num_genes, include_annotations):

    num_slim_genes = len(slim_mrcas_by_gene)
    if num_slim_genes == 0:
        avg_simulated_slim_Tc = 0
    else:
        avg_simulated_slim_Tc = sum(slim_mrcas_by_gene) / num_slim_genes

    if specks_mrcas_by_gene and (len(specks_mrcas_by_gene) > 0):
        num_specks_genes = len(specks_mrcas_by_gene)
        avg_simulated_specKS_Tc = sum(specks_mrcas_by_gene) / num_specks_genes
    else:
        avg_simulated_specKS_Tc = 0

    avg_simulated_specKS_Tc_in_years = avg_simulated_specKS_Tc * 10 ** 6
    if not xmax:
        xmax = max(slim_mrcas_by_gene)
    bins = np.arange(0, xmax, bin_size)

    two_Ne = 2.0 * Ne
    print("Tc plot bin_size_in_time: " + str(bin_size))
    kingman = [min(total_num_genes,
                   (bin_size * total_num_genes / two_Ne) * math.e ** ((-1 * i) / two_Ne))
               for i in bins]

    if slim_mrcas_by_gene:

        if include_annotations:
            my_label = 'DemographiKS Tcoal by gene\n' \
                       + "(" + str(num_slim_genes) + " genes in genome,\n" \
                       + "avg Tc " + str(int(avg_simulated_slim_Tc)) + " generations)"
        else:
            my_label = 'DemographiKS Tcoal by gene'

        bin_centers=[0.5*(bins[k]+bins[k+1]) for k in range(0,len(bins)-1)]

        n, my_bins, patches =this_ax.hist(slim_mrcas_by_gene, bins=bins, facecolor='b', alpha=0.25,
                     density=False, label=my_label)

        #this_ax.plot(bin_centers,n, color='b', alpha=0.25,
        #             label=my_label)

    if specks_mrcas_by_gene:
        num_specks_genes = len(specks_mrcas_by_gene)
        # note specks results are in millions of years,
        # so we have to convert millions of years to just years.
        million_years_to_years = 10 ** 6
        if include_annotations:
            my_label = 'SpecKS Tcoal by gene\n' \
                       + "(" + str(num_specks_genes) + " genes in genome),\n" \
                       + "avg Tc " + str(int(avg_simulated_specKS_Tc_in_years)) + " generations)"
        else:
            my_label = 'SpecKS Tcoal by gene'
        specks_mrcas_by_gene_in_YRs = [m * million_years_to_years for m in specks_mrcas_by_gene]

        n, my_bins, patches =this_ax.hist(specks_mrcas_by_gene_in_YRs, bins=bins,
                                          facecolor='c', alpha=0.25,
                     density=False, label=my_label)

        bin_centers=[0.5*(bins[k]+bins[k+1]) for k in range(0,len(bins)-1)]
        #this_ax.plot(bin_centers,n, color='b', alpha=0.25,
        #             label=my_label)

    #if theoretical_mrcas_by_gene:
    #    num_theory_genes = len(theoretical_mrcas_by_gene)
    #    avg_theory_Tc = sum(theoretical_mrcas_by_gene) / num_theory_genes
    #    this_ax.hist(theoretical_mrcas_by_gene, bins=bins, facecolor='r', alpha=0.25,
    #                 label='Theoretical Tcoal by gene\n'
    #                       + "(" + str(num_theory_genes) + " genes in genome,\n"
    #                       + "avg Tc " + str(int(avg_theory_Tc)) + " generations)",
    #                 density=False)

    this_ax.plot(bins, kingman, c='red', label='Expectations under Kingman,\n'
                                               + "exp. avg Tc " + str(int(two_Ne)) + " generations)")

    if ymax:
        this_ax.set(ylim=[0, ymax])
    this_ax.set(xlim=[0, xmax])
    this_ax.set(xlabel="MRCA time")
    this_ax.set(title=title)
    # this_ax.legend(loc='upper center')
    #this_ax.legend()
    return avg_simulated_slim_Tc
if __name__ == '__main__':
    unittest.main()

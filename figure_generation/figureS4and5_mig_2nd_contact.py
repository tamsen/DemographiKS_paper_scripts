import os
import unittest
import math
import sklearn
from matplotlib import pyplot as plt

from figure_generation.multi_row_Ks_aggregation_plots import make_multi_row_Ks_fig_with_subplots, SPKSvsDGKSPlotData


class MyTestCase(unittest.TestCase):

    def test_figS4_migration_as_secondary_contact(self):

        output_png_path='/home/tamsen/Data/DemographiKS_output_from_mesx/figS4_migration_as_secondary_contact_v2.png'
        qvalues_by_row=[0,1,1,1,1]
        r1_demographiKS_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Mg2_50Y'
        r1_run_list_name = " "
        r1_demographics_run_list = [
            'DGKS_0p00Mig50Y_figS4_row1.v3_m04d22y2026_h09m51s43',
            'DGKS_0p01Mig50Y_figS4_row1.v3_m04d22y2026_h09m53s05',
            'DGKS_0p1Mig50Y_figS4_row1.v3_m04d22y2026_h09m53s05',
            'DGKS_1p0Mig50Y_figS4_row1.v3_m04d22y2026_h09m53s05',
            'DGKS_10p0Mig50Y_figS4_row1.v3_m04d22y2026_h09m53s06',
        ]


        r1_specks_run_list = [False, False, False, False, False, False]
        r1_xmax_Ks = [0.02 for f in r1_demographics_run_list]
        r1_bins = [xmax_KS_i / 25 for xmax_KS_i in r1_xmax_Ks]
        r1_ymax_Ks = [250 for f in r1_demographics_run_list]
        r1_show_Ks_predictions = [False, False, False]
        r1_include_annotation = False
        r1_which_plot_panels_to_show_legend = [4]
        #r1_plot_title_lamda = lambda config: "Ks at Tnow\n" + \
        #                                  "Mig duration:" + \
        #                                  str(qvalues_by_row[1]*(config.mig_stop - config.mig_start)) + " gen"

        r1_plot_title_lamda = lambda config: "Ks at Tnow\n" + \
                                          "Mig duration: " + \
                                          str(qvalues_by_row[1]*(config.mig_stop - config.mig_start)) + " gen, " +\
                                          "rate: " +  str(qvalues_by_row[1]*(config.mig_rate)*100) + "%"

        row1_data = SPKSvsDGKSPlotData(r1_demographiKS_data_path, r1_demographics_run_list, r1_run_list_name,
                                       r1_demographiKS_data_path, r1_specks_run_list,
                                       r1_bins, r1_xmax_Ks, r1_ymax_Ks,
                                       r1_show_Ks_predictions,
                                       r1_include_annotation, r1_which_plot_panels_to_show_legend, r1_plot_title_lamda)


        r2_demographiKS_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Mg2_10KY'
        r2_run_list_name = " "
        r2_demographics_run_list = [
            'DGKS_0p00Mig10KY_figS5_row1.v3_m04d22y2026_h09m51s55',
            'DGKS_0p01Mig10KY_figS5_row1.v3_m04d22y2026_h09m52s55',
            'DGKS_0p1Mig10KY_figS5_row1.v3_m04d22y2026_h09m52s55',
            'DGKS_1p0Mig10KY_figS5_row1.v3_m04d22y2026_h09m52s55',
            'DGKS_10p0Mig10KY_figS5_row1.v3_m04d22y2026_h09m52s55',
            ]

        r2_plot_title_lamda = lambda config: "Ks at Tnow\n" + \
                                          "Mig duration: " + \
                                          str(qvalues_by_row[1]*(config.mig_stop - config.mig_start)) + " gen, " +\
                                          "rate: " +  str(qvalues_by_row[1]*(config.mig_rate)*100) + "%"


        row2_data = SPKSvsDGKSPlotData(r2_demographiKS_data_path, r2_demographics_run_list, r2_run_list_name,
                                       r1_demographiKS_data_path, r1_specks_run_list,
                                       r1_bins, r1_xmax_Ks, r1_ymax_Ks,
                                       r1_show_Ks_predictions,
                                       r1_include_annotation, r1_which_plot_panels_to_show_legend, r2_plot_title_lamda)


        make_multi_row_Ks_fig_with_subplots([row1_data,row2_data],
                                            output_png_path)




    def test_Ks_diffs_for_50_gen_of_Migration(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Mg2_50Y'
        data_file = "Ks_for_50yr_of_Mig_with_varying_rates_fig_Fig_S4.csv"
        png_out = "Ks_for_50yr_of_Mig_with_varying_rates_fig_Fig_S4_overlay.png"
        rmse_plot_label = "Ks perturbation vs mig rate (50 years contact)"
        make_RMSE_mig_plot(data_file, demographiKS_out_path, png_out, rmse_plot_label)

    def test_Ks_diffs_for_10K_gen_of_Migration(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Mg2_10KY'
        data_file = "Ks_for_10Kyr_of_Mig_with_varying_ratesfig_Fig_S5.csv"
        png_out = "Ks_for_10Kyr_of_Mig_with_varying_ratesfig_Fig_S5_overlay.png"
        rmse_plot_label = "Ks perturbation vs mig rate (10000 years contact)"
        make_RMSE_mig_plot(data_file, demographiKS_out_path, png_out, rmse_plot_label)


def make_RMSE_mig_plot(data_file, demographiKS_out_path, png_out, rmse_plot_label):
    with open(os.path.join(demographiKS_out_path, data_file), 'r') as f:
        lines = f.readlines()
    no_mig_dat_splat = lines[0].split(",")
    bins = no_mig_dat_splat[1].split(" ")
    bin_floats = [float(b) for b in bins]
    xs = [0.5 * (bin_floats[i] + bin_floats[i + 1]) for i in range(0, len(bins) - 1)]
    ys_no_mig = [float(d) for d in no_mig_dat_splat[2].split(" ")]
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    percent_styles = ["o", "x", "+", ">", "s"]
    rmses = []
    mig_rates = []
    p_int=0
    for line in lines:
        line_splat = line.split(",")
        line_label = line_splat[0]
        ys = [float(d) for d in line_splat[2].split(" ")]

        diffs = [ys[j] - ys_no_mig[j] for j in range(0, len(ys))]
        rmse = math.sqrt(sum([d * d for d in diffs]) / len(diffs))

        ax[0].plot(xs, ys, label=line_label)
        rmses.append(rmse)
        if "False" in line_label:
            mig_rate = 0
        else:
            mig_rate = float(line_label.replace("Mig rate:", ""))

        ax[1].scatter(mig_rate, rmse, label=line_label, marker = percent_styles[p_int])
        mig_rates.append(mig_rate)
        p_int =p_int+1
    ax[1].plot(mig_rates, rmses, color='gray')
    ax[0].legend()
    ax[1].legend()
    # this_ax.set(xlim=[0, xmax])
    ax[0].set(xlabel="Ks")
    ax[0].set(ylabel="# paralogs")
    ax[1].set(xlabel="Mig rate")
    ax[1].set(ylabel="Ks perturbation (RMSE)")
    ax[0].set(title="Ks histogram overlay")
    ax[1].set(title=rmse_plot_label)
    plt.savefig(os.path.join(demographiKS_out_path, png_out))
    plt.clf()


if __name__ == '__main__':
    unittest.main()

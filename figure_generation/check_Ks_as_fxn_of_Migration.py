import math
import os
import unittest

from matplotlib import pyplot as plt

from figure_generation.colors_for_figures import lighten_color
from figure_generation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots


class TestKsByMig(unittest.TestCase):

    def test_Ks_for_5_gen_of_Migration(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_Mg/5years'

        run_list_name = "Ks_for_5yr_of_Mig_with_varying_rates_fig_Fig_R-M5."
        demographics_run_list = [False,
                                'Mig14v4_m02d14y2025_h09m46s25','Mig15v4_m02d15y2025_h17m50s44',
                                 'Mig16v4_m02d15y2025_h17m50s37','Mig17v4_m02d15y2025_h17m50s41']

        specks_TE5_run_list = [False, False, False, False, False,False,False]
        xmax_Ks = [0.02 for f in demographics_run_list]
        bin_sizes_Ks = [xmax_KS_i/25 for xmax_KS_i in xmax_Ks]
        xmax_Tc = [10000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/25 for xmax_Tc_i in xmax_Tc]
        ymax_KS = [2000 for f in demographics_run_list]
        ymax_Tc = [2400 for f in demographics_run_list]

        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation=False
        plot_title_lamda = lambda config: "Ks at Tnow\n"+ "Mig rate:" + str(config.mig_rate)
        plots_to_show_legend=[1,2,3,4,5]
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, demographiKS_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions, include_annotation,
                                     plots_to_show_legend,
                                     plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here

    def test_Ks_for_100_gen_of_Migration(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_Mg/100years'

        run_list_name = "Ks_for_100yr_of_Mig_with_varying_ratesfig_Fig_R-M6."
        demographics_run_list = [False,
                                 'Mig14v4_m02d14y2025_h09m46s25',
                                'Mig13v4_m02d14y2025_h09m46s25',
                                 'Mig12v4_m02d14y2025_h09m46s25',
                                 'Mig11v4_m02d14y2025_h09m46s25']
        specks_TE5_run_list = [False, False, False, False, False, False, False]


        xmax_Ks = [0.02 for f in demographics_run_list]
        bin_sizes_Ks = [xmax_KS_i / 25 for xmax_KS_i in xmax_Ks]
        xmax_Tc = [10000 for f in demographics_run_list]
        bin_sizes_Tc = [xmax_Tc_i / 25 for xmax_Tc_i in xmax_Tc]
        ymax_KS = [2000 for f in demographics_run_list]
        ymax_Tc = [2400 for f in demographics_run_list]

        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation=False
        plots_to_show_legend=[1,2,3,4,5]
        plot_title_lamda = lambda config: "Ks at Tnow\n"+ "Mig rate:" + str(config.mig_rate)
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, demographiKS_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions,
                                     include_annotation,plots_to_show_legend,
                                     plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here

    def test_Ks_for_gradual_speciation_0p1_Migration(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_Mg/0p1percent'
        run_list_name = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_0p1percent"
        demographics_run_list = [False,
            'Mig14v4_m02d14y2025_h09m46s25',
            'Mig26v4p001_m02d19y2025_h09m33s24',
            'Mig27v4p001_m02d19y2025_h09m33s21',
            'Mig28v4p001_m02d19y2025_h09m33s19',
            'Mig29v4p001_m02d19y2025_h09m33s17']

        specks_TE5_run_list = [False, False, False, False, False,False,False]
        xmax_Ks = [0.02 for f in demographics_run_list]
        bin_sizes_Ks = [xmax_KS_i/25 for xmax_KS_i in xmax_Ks]
        xmax_Tc = [10000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/25 for xmax_Tc_i in xmax_Tc]
        ymax_KS = [2000 for f in demographics_run_list]
        ymax_Tc = [2400 for f in demographics_run_list]
        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation=False
        plots_to_show_legend=[]
        plot_title_lamda = lambda config: "Ks at Tnow\n"+ \
                                          "Mig duration:" + \
                                          str(config.mig_stop-config.mig_start) + " gen"
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, demographiKS_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions,
                                     include_annotation, plots_to_show_legend,
                                     plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here

    def test_Ks_for_gradual_speciation_1p0_Migration(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_Mg/1percent'
        run_list_name = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_1percent."
        demographics_run_list = [False,
         'Mig14v4_m02d14y2025_h09m46s25','Mig26v4p1_m02d17y2025_h18m42s39',
        'Mig27v4p1_m02d17y2025_h18m42s40',
         'Mig28v4p1_m02d17y2025_h18m42s43', 'Mig29v4p1_m02d17y2025_h18m42s46']

        specks_TE5_run_list = [False, False, False, False, False,False,False]
        xmax_Ks = [0.02 for f in demographics_run_list]
        bin_sizes_Ks = [xmax_KS_i/25 for xmax_KS_i in xmax_Ks]
        xmax_Tc = [10000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/25 for xmax_Tc_i in xmax_Tc]
        ymax_KS = [2000 for f in demographics_run_list]
        ymax_Tc = [2400 for f in demographics_run_list]
        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation=False
        plots_to_show_legend=[]
        plot_title_lamda = lambda config: "Ks at Tnow\n" + \
                                          "Mig duration:" + \
                                          str(config.mig_stop - config.mig_start)+ " gen"
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, demographiKS_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions,
                                     include_annotation, plots_to_show_legend,
                                     plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here


    def test_Ks_for_gradual_speciation_10_Migration(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_Mg/10percent'
        run_list_name = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_10percent."
        demographics_run_list = [False,
                                 'Mig14v4_m02d14y2025_h09m46s25',
                                 'Mig26v4_m02d16y2025_h14m03s20',
                                 'Mig27v4_m02d16y2025_h14m03s22',
                                 'Mig28v4_m02d16y2025_h14m03s34',
                                 'Mig29v4_m02d16y2025_h14m03s57']

        specks_TE5_run_list = [False, False, False, False, False,False,False]
        xmax_Ks = [0.02 for f in demographics_run_list]
        bin_sizes_Ks = [xmax_KS_i/25 for xmax_KS_i in xmax_Ks]
        xmax_Tc = [10000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/25 for xmax_Tc_i in xmax_Tc]
        ymax_KS = [2000 for f in demographics_run_list]
        ymax_Tc = [2400 for f in demographics_run_list]
        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation=False
        plots_to_show_legend=[]
        plot_title_lamda = lambda config: "Ks at Tnow\n"+ \
                                          "Mig duration:" + \
                                          str(config.mig_stop-config.mig_start)+ " gen"
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, demographiKS_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions,
                                     include_annotation,
                                     plots_to_show_legend, plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here

    def test_Ks_diffs_for_5_gen_of_Migration(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_Mg/5years'
        data_file = "Ks_for_5yr_of_Mig_with_varying_rates_fig_Fig_R-M5_save.csv"
        png_out = "Ks_for_5yr_of_Mig_with_varying_rates_fig_Fig_R-M5_overlay.png"
        rmse_plot_label = "Ks perturbation vs mig rate (5 years contact)"
        make_RMSE_mig_plot(data_file, demographiKS_out_path, png_out, rmse_plot_label)

    def test_Ks_diffs_for_100_gen_of_Migration(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_Mg/100years'
        data_file = "Ks_for_100yr_of_Mig_with_varying_ratesfig_Fig_R-M6_save.csv"
        png_out = "Ks_for_100yr_of_Mig_with_varying_ratesfig_Fig_R-M6_overlay.png"
        rmse_plot_label = "Ks perturbation vs mig rate (100 years contact)"
        make_RMSE_mig_plot(data_file, demographiKS_out_path, png_out, rmse_plot_label)

    def test_Ks_diffs_for_gradual_speciation(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_Mg'
        data_file_0p1 = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_0p1percent.csv"
        data_file_1p0 = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_1percent.csv"
        data_file_10 = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_10percent.csv"
        png_out = "Ks_for_gradual_speciation_overlay.png"
        rmse_plot_label = "Ks perturbation vs mig rate"

        with open(os.path.join(demographiKS_out_path, data_file_0p1), 'r') as f:
            lines_0p1 = f.readlines()
        with open(os.path.join(demographiKS_out_path, data_file_1p0), 'r') as f:
            lines_1p0 = f.readlines()
        with open(os.path.join(demographiKS_out_path, data_file_10), 'r') as f:
            lines_10 = f.readlines()

        data_list=[lines_0p1,lines_1p0,lines_10]

        no_mig_dat_splat = lines_0p1[0].split(",")
        bins = no_mig_dat_splat[1].split(" ")
        bin_floats = [float(b) for b in bins]
        xs = [0.5 * (bin_floats[i] + bin_floats[i + 1]) for i in range(0, len(bins) - 1)]
        ys_no_mig = [float(d) for d in no_mig_dat_splat[2].split(" ")]
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        #rmses = []
        mig_rates = ["0.1%","1%","10%"]
        color_adjustment=[1.0,.8,.6]
        marker_style = ["o","x","+",">","s"]
        #line_style = ["-", ":", "--","-.","-."]
        colors=['b','g','r']
        duration_colors=['b','g','r','c','purple']
        for percent_index in range(0,len(data_list)):
            lines_p=data_list[percent_index]
            percent_color=colors[percent_index]
            rmses = []
            mig_durations = []
            for duration_index in range(0,len(lines_p)):
                line=lines_p[duration_index]
                line_splat = line.split(",")
                duration_string=line_splat[0].replace("Mig duration:", "").replace(" gen", "")
                line_label = "mig rate:" +mig_rates[percent_index]
                ys = [float(d) for d in line_splat[2].split(" ")]

                diffs = [ys[j] - ys_no_mig[j] for j in range(0, len(ys))]
                rmse = math.sqrt(sum([d * d for d in diffs]) / len(diffs))

                #base_color=duration_colors[duration_index]
                #adjusted_color=lighten_color(base_color, amount=color_adjustment[percent_index])
                #if duration_index==0:
                #    adjusted_color =duration_colors[duration_index]
                #    ax[0].plot(xs, ys, label=line_label,color=adjusted_color,
                #               linestyle=line_style[duration_index])
                #else:
                #    ax[0].plot(xs, ys, color=adjusted_color,
                #               linestyle=line_style[duration_index])

                if duration_index==0:
                    ax[0].plot(xs, ys, label=line_label,color=percent_color,
                               marker=marker_style[duration_index],
                               mfc='gray', mec='gray', markersize=5,alpha=0.5)
                else:
                    ax[0].plot(xs, ys, color=percent_color,
                               marker=marker_style[duration_index],
                               mfc='gray', mec='gray', markersize=5,alpha=0.5)
                duration_float=float(duration_string)
                rmses.append(rmse)
                mig_durations.append(duration_float)
                if percent_index==0:
                    ax[1].scatter(duration_float, rmse,
                              marker=marker_style[duration_index],
                              color='gray', s=50, label=duration_string + " gen")

            ax[1].plot(mig_durations, rmses,color=percent_color,alpha=0.5)

        ax[0].legend()
        ax[1].legend()
        # this_ax.set(xlim=[0, xmax])
        ax[0].set(xlabel="Ks")
        ax[0].set(ylabel="# paralogs")
        ax[1].set(xlabel="Mig duration")
        ax[1].set(ylabel="Ks perturbation (RMSE)")
        ax[0].set(title="Ks histogram overlay")
        ax[1].set(title=rmse_plot_label)
        plt.savefig(os.path.join(demographiKS_out_path, png_out))
        plt.clf()

    def test_Ks_diffs_for_gradual_speciation_try_again(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_Mg'
        data_file_0p1 = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_0p1percent.csv"
        data_file_1p0 = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_1percent.csv"
        data_file_10 = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_10percent.csv"
        png_out = "Ks_for_gradual_speciation_overlay_alternative.png"
        rmse_plot_label = "Ks perturbation vs mig rate"
        marker_displacement=[0,0.1,-0.1,0.2,-0.2]
        with open(os.path.join(demographiKS_out_path, data_file_0p1), 'r') as f:
            lines_0p1 = f.readlines()
        with open(os.path.join(demographiKS_out_path, data_file_1p0), 'r') as f:
            lines_1p0 = f.readlines()
        with open(os.path.join(demographiKS_out_path, data_file_10), 'r') as f:
            lines_10 = f.readlines()

        data_list=[lines_0p1,lines_1p0,lines_10]

        no_mig_dat_splat = lines_0p1[0].split(",")
        bins = no_mig_dat_splat[1].split(" ")
        bin_floats = [float(b) for b in bins]
        xs = [0.5 * (bin_floats[i] + bin_floats[i + 1]) for i in range(0, len(bins) - 1)]
        ys_no_mig = [float(d) for d in no_mig_dat_splat[2].split(" ")]
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        #rmses = []
        mig_rates = ["0.1%","1%","10%"]
        color_adjustment=[1.0,.8,.6]
        percent_styles = ["o","x","+",">","s"]
        #line_style = ["-", ":", "--","-.","-."]
        #colors=['b','g','r']
        duration_colors=['b','g','r','c','purple']
        for percent_index in range(0,len(data_list)):
            lines_p=data_list[percent_index]
            #percent_color=colors[percent_index]
            rmses = []
            mig_durations = []
            for duration_index in range(0,len(lines_p)):
                line=lines_p[duration_index]
                line_splat = line.split(",")
                duration_string=line_splat[0].replace("Mig duration:", "").replace(" gen", "")
                percent_string = "mig rate:" +mig_rates[percent_index]
                ys = [float(d) for d in line_splat[2].split(" ")]

                diffs = [ys[j] - ys_no_mig[j] for j in range(0, len(ys))]
                rmse = math.sqrt(sum([d * d for d in diffs]) / len(diffs))

                if percent_index==0:
                    ax[0].plot(xs, ys, label=duration_string + " gen"
                               ,color=duration_colors[duration_index],
                               marker=percent_styles[percent_index],
                               mfc='gray', mec='gray', markersize=5)
                else:
                    ax[0].plot(xs, ys, color=duration_colors[duration_index],
                               marker=percent_styles[percent_index],
                               mfc='gray', mec='gray', markersize=5)
                duration_float=float(duration_string)
                rmses.append(rmse)
                mig_durations.append(duration_float)
                if duration_index==0:
                    ax[1].scatter(duration_float, rmse,
                              marker=percent_styles[percent_index],
                              color='gray',
                                  s=50, label=percent_string, zorder=0)
                ax[1].scatter(duration_float,rmse +500*marker_displacement[percent_index],
                              marker=percent_styles[percent_index],
                              color=duration_colors[duration_index],
                                  s=50)

            ax[1].plot(mig_durations,[r+500*marker_displacement[percent_index] for r in rmses],
                       color='gray',
                       alpha=0.5)

        ax[0].legend()
        ax[1].legend()
        # this_ax.set(xlim=[0, xmax])
        ax[0].set(xlabel="Ks")
        ax[0].set(ylabel="# paralogs")
        ax[1].set(xlabel="Mig duration")
        ax[1].set(ylabel="Ks perturbation (RMSE)")
        ax[0].set(title="Ks histogram overlay")
        ax[1].set(title=rmse_plot_label)
        plt.savefig(os.path.join(demographiKS_out_path, png_out))
        plt.clf()


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
        mig_rate=float(line_label.replace("Mig rate:", ""))
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

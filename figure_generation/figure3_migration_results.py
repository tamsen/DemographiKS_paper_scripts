import os
import unittest
import math
import sklearn
from matplotlib import pyplot as plt

from figure_generation.multi_row_Ks_aggregation_plots import make_multi_row_Ks_fig_with_subplots, SPKSvsDGKSPlotData


class MyTestCase(unittest.TestCase):

    def test_fig3_migration_as_gradual_speciation(self):

        output_png_path='/home/tamsen/Data/DemographiKS_output_from_mesx/fig3_migration_as_gradual_speciation.png'
        qvalues_by_row=[0,100,100,10,1]
        r1_demographiKS_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Mg2/0p01percent'
        r1_run_list_name = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_0p01percent_v2."
        r1_demographics_run_list = [
            'DGKS_0p0Mig0KY_fig2_row1.v1_m04d02y2026_h14m14s39',
                                 'DGKS_0p01Mig50KY_fig2_row1.v1_m04d02y2026_h18m14s58',
                                 'DGKS_0p01Mig100KY_fig2_row1.v1_m04d02y2026_h18m14s58',
                                 'DGKS_0p01Mig250KY_fig2_row1.v1_m04d02y2026_h18m14s58',
                                 'DGKS_0p01Mig500KY_fig2_row1.v1_m04d02y2026_h18m14s58',
                                 ]
        r1_specks_run_list = [False, False, False, False, False, False, False]
        r1_xmax_Ks = [0.02 for f in r1_demographics_run_list]
        r1_bins = [xmax_KS_i / 25 for xmax_KS_i in r1_xmax_Ks]
        r1_ymax_Ks = [250 for f in r1_demographics_run_list]
        r1_show_Ks_predictions = [False, False, False]
        r1_include_annotation = False
        r1_which_plot_panels_to_show_legend = []
        r1_plot_title_lamda = lambda config: "Ks at Tnow\n" + \
                                          "Mig duration:" + \
                                          str(qvalues_by_row[1]*(config.mig_stop - config.mig_start)) + " gen"

        row1_data = SPKSvsDGKSPlotData(r1_demographiKS_data_path, r1_demographics_run_list, r1_run_list_name,
                                       r1_demographiKS_data_path, r1_specks_run_list,
                                       r1_bins, r1_xmax_Ks, r1_ymax_Ks,
                                       r1_show_Ks_predictions,
                                       r1_include_annotation, r1_which_plot_panels_to_show_legend, r1_plot_title_lamda)

        r2_demographiKS_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Mg2/0p10percent'
        r2_run_list_name = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_0p1percent.v2."
        r2_demographics_run_list = ['DGKS_0p0Mig0KY_fig2_row1.v1_m04d02y2026_h14m14s39',
                                     'DGKS_0p1Mig50KY_fig2_row1.v1_m04d02y2026_h14m12s46',
                                     'DGKS_0p1Mig100KY_fig2_row1.v1_m04d02y2026_h14m12s46',
                                     'DGKS_0p1Mig250KY_fig2_row1.v1_m04d02y2026_h14m12s46',
                                     'DGKS_0p1Mig500KY_fig2_row1.v1_m04d02y2026_h14m12s47']
        r2_plot_title_lamda = lambda config: "Ks at Tnow\n" + \
                                          "Mig duration:" + \
                                          str(qvalues_by_row[2]*(config.mig_stop - config.mig_start)) + " gen"

        row2_data = SPKSvsDGKSPlotData(r2_demographiKS_data_path, r2_demographics_run_list, r2_run_list_name,
                                       r1_demographiKS_data_path, r1_specks_run_list,
                                       r1_bins, r1_xmax_Ks, r1_ymax_Ks,
                                       r1_show_Ks_predictions,
                                       r1_include_annotation, r1_which_plot_panels_to_show_legend, r2_plot_title_lamda)

        r3_demographiKS_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Mg2/1p00percent'
        r3_run_list_name = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_1percent_v2."
        r3_demographics_run_list = ['DGKS_0p0Mig0KY_fig2_row1.v1_m04d02y2026_h14m14s39',
                                     'DGKS_1p0Mig50KY_fig2_row1.v1_m04d02y2026_h14m44s20',
                                     'DGKS_1p0Mig100KY_fig2_row1.v1_m04d02y2026_h14m48s10',
                                     'DGKS_1p0Mig250KY_fig2_row1.v1_m04d02y2026_h14m48s14',
                                     'DGKS_1p0Mig500KY_fig2_row1.v1_m04d02y2026_h14m48s19',
                                     ]

        r3_plot_title_lamda = lambda config: "Ks at Tnow\n" + \
                                          "Mig duration:" + \
                                          str(qvalues_by_row[3]*(config.mig_stop - config.mig_start)) + " gen"

        row3_data = SPKSvsDGKSPlotData(r3_demographiKS_data_path, r3_demographics_run_list, r3_run_list_name,
                                       r1_demographiKS_data_path, r1_specks_run_list,
                                       r1_bins, r1_xmax_Ks, r1_ymax_Ks,
                                       r1_show_Ks_predictions,
                                       r1_include_annotation, r1_which_plot_panels_to_show_legend, r3_plot_title_lamda)


        r4_demographiKS_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Mg2/10percent'
        r4_run_list_name = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_10percent_v2."
        r4_demographics_run_list = ['DGKS_0p0Mig0KY_fig2_row1.v1_m04d02y2026_h14m14s39',
                                     'DGKS_10p0Mig50KY_fig2_row1.v1_m04d02y2026_h18m08s27',
                                     'DGKS_10p0Mig100KY_fig2_row1.v1_m04d02y2026_h18m08s27',
                                     'DGKS_10p0Mig250KY_fig2_row1.v1_m04d02y2026_h18m08s27',
                                     'DGKS_10p0Mig500KY_fig2_row1.v1_m04d02y2026_h18m08s27',
                                     ]

        r4_plot_title_lamda = lambda config: "Ks at Tnow\n" + \
                                          "Mig duration:" + \
                                          str(qvalues_by_row[4]*(config.mig_stop - config.mig_start)) + " gen"
        row4_data = SPKSvsDGKSPlotData(r4_demographiKS_data_path, r4_demographics_run_list, r4_run_list_name,
                                       r1_demographiKS_data_path, r1_specks_run_list,
                                       r1_bins, r1_xmax_Ks, r1_ymax_Ks,
                                       r1_show_Ks_predictions,
                                       r1_include_annotation, r1_which_plot_panels_to_show_legend, r4_plot_title_lamda)

        make_multi_row_Ks_fig_with_subplots([row1_data,row2_data,row3_data,row4_data],
                                            output_png_path)


    def test_Ks_diffs_for_gradual_speciation_nicer_plot(self):

        output_png_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/fig3_migration_overlay.png'
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Ks_vs_Mg2'
        data_file_0p01 = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_0p01percent_v2._out.csv"
        data_file_0p1 = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_0p1percent_v2._out.csv"
        data_file_1p0 = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_1percent_v2._out.csv"
        data_file_10 = "Ks_for_gradual_speciation_ie_Mig_with_varying_duration_10percent_v2._out.csv"
        png_out = "Ks_for_gradual_speciation_overlay_alternative_v2.png"
        #rmse_plot_label = "Ks perturbation vs mig rate"
        #marker_displacement=[0,0.1,-0.1,0.2,-0.2]
        marker_displacement = [0, 0, 0,0, 0, 0]
        qvalues_by_row = [0, 100, 100, 10, 1]

        with open(os.path.join(demographiKS_out_path,"data_files", data_file_0p01), 'r') as f:
            lines_0p01 = f.readlines()
        with open(os.path.join(demographiKS_out_path, "data_files",data_file_0p1), 'r') as f:
            lines_0p1 = f.readlines()
        with open(os.path.join(demographiKS_out_path, "data_files",data_file_1p0), 'r') as f:
            lines_1p0 = f.readlines()
        with open(os.path.join(demographiKS_out_path, "data_files",data_file_10), 'r') as f:
            lines_10 = f.readlines()

        data_list=[lines_0p01,lines_0p1,lines_1p0,lines_10]

        #data_qscaling  no longer needed, this is taken care of by the lamda
        data_qscaling = [1, 1, 1, 1]

        no_mig_dat_splat = lines_0p1[0].split(",")
        bins = no_mig_dat_splat[1].split(" ")
        bin_floats = [float(b) for b in bins]
        xs = bin_floats[0:-1]
        #xs = [0.5 * (bin_floats[i] + bin_floats[i + 1]) for i in range(0, len(bins) - 1)]
        ys_no_mig = [float(d) for d in no_mig_dat_splat[2].split(" ")]
        fig, ax = plt.subplots(1, 2, figsize=(12, 4), dpi=500)

        mig_rates = ["0.01%","0.1%","1%","10%"]
        color_adjustment=[1.0,.8,.6]
        percent_styles = ["o","x","+",">","s","<"]
        #line_style = ["-", ":", "--","-.","-."]
        #colors=['b','g','r']
        duration_colors=['b','g','r','c','purple','y']
        percent_alphas=[0.1,0.2,0.5,1.0]
        for percent_index in range(0,len(data_list)):
            lines_p=data_list[percent_index]
            data_qscaling_for_percent=data_qscaling[percent_index]
            #percent_color=colors[percent_index]
            rmses = []
            mig_durations = []
            r2s=[]
            for duration_index in range(0,len(lines_p)):
                line=lines_p[duration_index]
                line_splat = line.split(",")
                duration_string=line_splat[0].replace("Mig duration:", "").replace(" gen", "")
                percent_string = "mig rate:" +mig_rates[percent_index]
                ys = [float(d) for d in line_splat[2].split(" ")]
                print("duration string: " + str(duration_string))
                diffs = [ys[j] - ys_no_mig[j] for j in range(0, len(ys))]
                rmse = math.sqrt(sum([d * d for d in diffs]) / len(diffs))
                if duration_index == 0:
                    ax1label =percent_string
                else:
                    ax1label = None

                if percent_index==0:
                    print(str(duration_colors))
                    print(str(duration_index))
                    ax[0].plot(xs, ys, label=duration_string + " gen"
                               ,color=duration_colors[duration_index])
                    #           marker=percent_styles[percent_index],
                    #           mfc='gray', mec='gray', markersize=5)

                    #ax[1].plot(xs, ys, label=ax1label
                    #           ,color=duration_colors[duration_index],
                    #           alpha=percent_alphas[percent_index])
                    #marker=percent_styles[percent_index],
                    #           mfc='gray', mec='gray', markersize=5,

                else:
                    ax[0].plot(xs, ys, color=duration_colors[duration_index],

                               alpha=percent_alphas[percent_index])
                    #ax[1].plot(xs, ys, label=ax1label
                    #           ,color=duration_colors[duration_index],
                    #           )
                    #           marker=percent_styles[percent_index],
                    #                                mfc='gray', mec='gray', markersize=5,

                duration_float = float(duration_string) * data_qscaling_for_percent
                print(" -- + --")
                print("raw mig duration:" +str(line_splat[0]))
                print("data_qscaling_for_percent: " + str( data_qscaling_for_percent))
                print("duration_string: " + str(duration_string))
                print("duration_float: " + str(duration_float))
                print(" -- - --")
                #r2_score = sklearn.metrics.r2_score(ys,ys_no_mig)
                #r2s.append(r2_score)
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
        #ax[1].set(ylabel="Ks perturbation (RMSE)")
        ax[0].set(title="Ks overlay by migration rate")
        ax[1].set(title="Ks perturbation due to migration, as RMSE")
        #ax[1].set(title=rmse_plot_label)
        plt.savefig(output_png_path)
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

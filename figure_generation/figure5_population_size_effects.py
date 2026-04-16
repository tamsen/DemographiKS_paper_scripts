import os
import unittest

from figure_generation.multirow_aggregations_auto_vs_allo import AUTOvsALLOPlotData, \
    make_multirow_Allo_vs_Auto_fig_with_subplots


class MyTestCase(unittest.TestCase):

    def test_fig5_population_size_effects(self):

        output_root_path = '/home/tamsen/Data/DemographiKS_output_from_mesx'
        output_png_path = os.path.join(output_root_path, 'fig5_population_size_effects.png')
        qvalues_by_row = [0, 1000, 1000]
        r1_auto_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_v2/NbVaries/NaFixedAt10K_Q1000'
        r1_allo_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_v2/NbVaries/NaFixedAt10K_Q1000'

        r1_allo_run_list = [
            'Allo_10KNa_10Nb_v4_m04d08y2026_h12m26s56',
            'Allo_10KNa_50Nb_v4_m04d08y2026_h12m26s56',
            'Allo_10KNa_100Nb_v4_m04d08y2026_h12m26s56',
            'Allo_10KNa_500Nb_v4_m04d08y2026_h12m26s58',
            'Allo_10KNa_500Nb_old_v4_m04d08y2026_h16m17s04']

        r1_auto_run_list = [
            'Auto_10KNa_10Nb_v4_m04d08y2026_h12m29s29',
            'Auto_10KNa_50Nb_v4_m04d08y2026_h12m29s29',
            'Auto_10KNa_100Nb_v4_m04d08y2026_h12m29s29',
            'Auto_10KNa_500Nb_v4_m04d08y2026_h12m29s29',
            'Auto_10KNa_500Nb_old_v4_m04d08y2026_h16m16s59']

        r1_xmax_Ks_array = [0.02 for f in r1_allo_run_list]
        r1_bin_sizes_Ks_array = [xmax_KS_i/25 for xmax_KS_i in r1_xmax_Ks_array]
        r1_ymax_Ks_array = [200 for f in r1_allo_run_list]
        r1_auto_run_list_name = None#"Ks_for_Allo_and_Auto_varying_varying_Nb_10K_Na_Fig_R-NaNb6_v2"
        r1_show_Ks_predictions = [False, False, False]
        r1_suptitle = "Auto and Allo Ks histograms\n"
        r1_plot_title_lamda = lambda config: (
                "Na = " + str(qvalues_by_row[1]*config.ancestral_Ne) +
                ", Nb = " + str(qvalues_by_row[1] * config.bottleneck_Ne) +
                "\nTwgd = " + str(int(qvalues_by_row[1]*config.WGD_time_Ge)))
        r1_which_plot_panels_to_show_legend = []
        r1_include_annotation = False

        row1_data = AUTOvsALLOPlotData(r1_auto_data_path, r1_auto_run_list, r1_auto_run_list_name,
                                       r1_allo_data_path, r1_allo_run_list,
                                       r1_bin_sizes_Ks_array, r1_xmax_Ks_array, r1_ymax_Ks_array,
                                       r1_show_Ks_predictions,
                                       r1_include_annotation,
                                       r1_which_plot_panels_to_show_legend, r1_plot_title_lamda)


        r2_auto_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_v2/NaVaries/NbFixedQ1000'
        r2_allo_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_v2/NaVaries/NbFixedQ1000'

        r2_allo_run_list = [
            'Allo_10KNb_10Na_v4_m04d08y2026_h17m18s47',
            'Allo_10KNb_50Na_v4_m04d08y2026_h17m18s47',
            'Allo_10KNb_100Na_v4_m04d08y2026_h17m18s47',
            'Allo_10KNb_500Na_v4_m04d09y2026_h10m44s19',
            'Allo_10KNa_500Nb_old_v4_m04d09y2026_h10m44s42']

        r2_auto_run_list = [
            'Auto_10KNa_10Nb_v4_m04d08y2026_h17m10s00',
            'Auto_10KNa_50Nb_v4_m04d08y2026_h17m10s00',
            'Auto_10KNa_100Nb_v4_m04d08y2026_h17m10s00',
            'Auto_10KNb_500Na_v4_m04d09y2026_h10m44s21',
            'Auto_10KNa_500Nb_old_v4_m04d09y2026_h10m44s51', ]

        r2_xmax_Ks_array = [0.02 for f in r2_auto_run_list]
        r2_bin_sizes_Ks_array = [xmax_KS_i / 25 for xmax_KS_i in r2_xmax_Ks_array]
        r2_ymax_Ks_array = [200 for f in r2_auto_run_list]
        r2_auto_run_list_name = None#"Ks_for_Allo_and_Auto_varying_varying_Na_10K_Nb_Fig_R-NaNb8_v2"
        r2_show_Ks_predictions = [False, False, False]
        r2_suptitle = "Auto and Allo Ks histograms\n"
        r2_include_annotation = False
        r2_which_plot_panels_to_show_legend=[]
        r2_plot_title_lamda = lambda config: (
                "Na = " + str(qvalues_by_row[1]*config.ancestral_Ne) +
                ", Nb = " + str(qvalues_by_row[1] * config.bottleneck_Ne) +
                "\nTwgd = " + str(int(qvalues_by_row[1]*config.WGD_time_Ge)))

        row2_data = AUTOvsALLOPlotData(r2_auto_data_path, r2_auto_run_list, r2_auto_run_list_name,
                                       r2_allo_data_path, r2_allo_run_list,
                                       r2_bin_sizes_Ks_array, r2_xmax_Ks_array, r2_ymax_Ks_array,
                                       r2_show_Ks_predictions,
                                       r2_include_annotation,
                                       r2_which_plot_panels_to_show_legend, r2_plot_title_lamda)

        r3_allo_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_v2/NaVaries/NaVariesNbfixed_OldWGD'
        r3_auto_data_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_v2/NaVaries/NaVariesNbfixed_OldWGD'

        r3_allo_run_list = [
            'Allo_10KNb_10Na_v5_m04d09y2026_h12m04s09',
            'Allo_10KNb_50Na_v5_m04d09y2026_h12m01s33',
            'Allo_10KNb_100Na_v5_m04d09y2026_h12m01s33',
            'Allo_10KNb_500Na_v5_m04d09y2026_h12m01s33',
            'Allo_10KNb_1000Na_v5_m04d09y2026_h12m01s35',
        ]

        r3_auto_run_list = [
            'Auto_10KNa_10Nb_v5_m04d09y2026_h12m10s10',
            'Auto_10KNa_50Nb_v5_m04d09y2026_h12m10s07',
            'Auto_10KNa_100Nb_v5_m04d09y2026_h12m09s56',
            'Auto_10KNb_500Na_v5_m04d09y2026_h12m09s59',
            'Auto_10KNb_1000Na_v5_m04d09y2026_h12m13s02',
        ]

        r3_xmax_Ks_array = [0.02 for f in r3_allo_run_list]
        r3_bin_sizes_Ks_array = [xmax_KS_i / 25 for xmax_KS_i in r3_xmax_Ks_array]

        r3_ymax_Ks_array = [200 for f in r3_allo_run_list]
        r3_auto_run_list_name = None #"Ks_for_Allo_and_Auto_varying_varying_Na_10K_Nb_Fig_R-NaNb8_olderWGD_v2"
        r3_show_Ks_predictions = [False, False, False]
        r3_suptitle = "Auto and Allo Ks histograms\n"
        r3_include_annotation = False
        r3_which_plot_panels_to_show_legend=[4]

        r3_plot_title_lamda = lambda config: (
                "Na = " + str(qvalues_by_row[1]*config.ancestral_Ne) +
                ", Nb = " + str(qvalues_by_row[1] * config.bottleneck_Ne) +
                "\nTwgd = " + str(int(qvalues_by_row[1]*config.WGD_time_Ge)))

        row3_data = AUTOvsALLOPlotData(r3_auto_data_path, r3_auto_run_list, r3_auto_run_list_name,
                                       r3_allo_data_path, r3_allo_run_list,
                                       r3_bin_sizes_Ks_array, r3_xmax_Ks_array, r3_ymax_Ks_array,
                                       r3_show_Ks_predictions,
                                       r3_include_annotation,
                                       r3_which_plot_panels_to_show_legend, r3_plot_title_lamda)

        make_multirow_Allo_vs_Auto_fig_with_subplots([row1_data, row2_data, row3_data], output_png_path)

if __name__ == '__main__':
    unittest.main()

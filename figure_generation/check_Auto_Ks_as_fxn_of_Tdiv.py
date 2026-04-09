import math
import os
import unittest

from matplotlib import pyplot as plt
from scipy.stats import expon

from figure_generation.ks_plot_aggregations_auto_vs_allo import make_Tc_Ks_Allo_vs_Auto_fig_with_subplots


class TestAlloVsAuto_Tdiv(unittest.TestCase):



    def test_allo_vs_auto_Ks_for_varying_varying_Tdiv_times_v2(self):


        demographiKS_auto_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_v2/TwgdVaries'
        demographiKS_allo_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_v2/TwgdVaries'

        demographics_allo_run_list=[False,
                                    'Allo_100KNa_100KMYA_v1_m04d07y2026_h10m45s12',
                                    'Allo_100KNa_200KMYA_v1_m04d07y2026_h09m05s10',
                                    'Allo_100KNa_500KMYA_v1_m04d07y2026_h08m53s38',
                                    'Allo_100KNa_1000KMYA_v1_m04d07y2026_h08m53s47']
        #'Auto_100KNa_100KMYA_v1_m04d06y2026_h17m53s29'
        demographics_auto_run_list=[False,
                                    'Auto_100KNa_100KMYA_v1_m04d06y2026_h17m53s29',
                                    'Auto_100KNa_200KMYA_v1_m04d07y2026_h09m05s43',
                                    'Auto_100KNa_500KMYA_v1_m04d06y2026_h17m58s12',
                                    'Auto_100KNa_1000KMYA_v1_m04d06y2026_h18m00s04']



        #xmax_Ks_allo = [0.05,0.025,0.025,0.025,0.025] #0.001  # max(demographiKS_ks_results)
        xmax_Ks_allo = [0.05, 0.05, 0.05, 0.05, 0.05]  # 0.001  # max(demographiKS_ks_results)
        xmax_Ks_auto = [0.02,0.02,0.02,0.02,0.02] #0.001  # max(demographiKS_ks_results)
        xmax_Tc = [20000 for f in demographics_allo_run_list]
        bin_sizes_Ks_allo = [xmax_Ks_i/25.0 for xmax_Ks_i in xmax_Ks_allo]
        bin_sizes_Ks_auto = [xmax_Ks_i /25.0 for xmax_Ks_i in xmax_Ks_auto]
        bin_sizes_Tc =  [xmax_Tc_i /25.0 for xmax_Tc_i in xmax_Tc]

        run_list_name="Ks_for_Allo_and_Auto_varying_varying_Tdiv_R-HomEx-Tdiv_v2"
        ymax_KS = [800,800,800,800,800]
        ymax_Tc = [400 for f in demographics_allo_run_list]
        #show_KS_predictions=[True,True,True]
        show_KS_predictions = [False, False, False]
        suptitle = "Allo and Auto Ks histograms\n"
        #                          "Recombination rate = 8e-9, Ne and BI constant"


        bin_sizes_Ks_array=[bin_sizes_Ks_auto,bin_sizes_Ks_allo]
        xmax_Ks_array=[xmax_Ks_auto,xmax_Ks_allo]
        ymax_Ks_array=[ymax_KS,ymax_KS]

        num_plot_rows =2
        include_annotation = False
        plot_title_lamda = lambda config: "Twgd&Tdiv:" + str(config.WGD_time_Ge) +" gen"
        make_Tc_Ks_Allo_vs_Auto_fig_with_subplots(num_plot_rows,bin_sizes_Ks_array, bin_sizes_Tc,
                                          demographiKS_allo_out_path,
                                                  demographics_allo_run_list, run_list_name,
                                          demographics_auto_run_list, demographiKS_auto_out_path,
                                     xmax_Ks_array, xmax_Tc, ymax_Ks_array,  ymax_Tc,suptitle,
                                     show_KS_predictions,include_annotation,
                                                  plot_title_lamda             )

        self.assertEqual(True, True)  # add assertion here


    def test_allo_vs_auto_Ks_for_varying_varying_Tdiv_times(self):


        demographiKS_auto_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto/Auto_vs_Tdiv'
        demographiKS_allo_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto/Auto_vs_Tdiv'
        #demographiKS_allo_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Tdiv'


        #specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Tdiv'
        demographics_allo_run_list=[False,
          'TE05fix__m01d06y2025_h11m17s15','TE07_fix__m01d08y2025_h15m09s22',
           'TE08_fix_m01d14y2025_h09m17s20','TE09_fix_m01d13y2025_h14m13s12']

        demographics_auto_run_list=[False,
                                    'Auto_Twgd10000_m01d29y2025_h17m22s23',
                                    'Auto_Twgd100000_m01d29y2025_h17m25s36',
                                    'Auto_Twgd500000_m01d29y2025_h17m22s23',
                                    'Auto_Twgd100000_m01d29y2025_h17m25s36']



        #xmax_Ks_allo = [0.05,0.025,0.025,0.025,0.025] #0.001  # max(demographiKS_ks_results)
        xmax_Ks_allo = [0.05, 0.05, 0.05, 0.05, 0.05]  # 0.001  # max(demographiKS_ks_results)
        xmax_Ks_auto = [0.02,0.02,0.02,0.02,0.02] #0.001  # max(demographiKS_ks_results)
        xmax_Tc = [12000 for f in demographics_allo_run_list]
        bin_sizes_Ks_allo = [xmax_Ks_i/25.0 for xmax_Ks_i in xmax_Ks_allo]
        bin_sizes_Ks_auto = [xmax_Ks_i /25.0 for xmax_Ks_i in xmax_Ks_auto]
        bin_sizes_Tc =  [xmax_Tc_i /25.0 for xmax_Tc_i in xmax_Tc]

        run_list_name="Ks_for_Allo_and_Auto_varying_varying_Tdiv_R-HomEx-Tdiv_test"
        ymax_KS = [800,800,800,800,800]
        ymax_Tc = [400 for f in demographics_allo_run_list]
        #show_KS_predictions=[True,True,True]
        show_KS_predictions = [False, False, False]
        suptitle = "Allo and Auto Ks histograms\n"
        #                          "Recombination rate = 8e-9, Ne and BI constant"


        bin_sizes_Ks_array=[bin_sizes_Ks_auto,bin_sizes_Ks_allo]
        xmax_Ks_array=[xmax_Ks_auto,xmax_Ks_allo]
        ymax_Ks_array=[ymax_KS,ymax_KS]

        num_plot_rows =2
        include_annotation = False
        plot_title_lamda = lambda config: "Twgd&Tdiv:" + str(config.WGD_time_Ge) +" gen"
        make_Tc_Ks_Allo_vs_Auto_fig_with_subplots(num_plot_rows,bin_sizes_Ks_array, bin_sizes_Tc,
                                          demographiKS_allo_out_path,
                                                  demographics_allo_run_list, run_list_name,
                                          demographics_auto_run_list, demographiKS_auto_out_path,
                                     xmax_Ks_array, xmax_Tc, ymax_Ks_array,  ymax_Tc,suptitle,
                                     show_KS_predictions,include_annotation,
                                                  plot_title_lamda             )

        self.assertEqual(True, True)  # add assertion here

    def test_expected_vs_sim_for_auto_vs_allo_Tdiv(self):

        output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/Auto/Auto_vs_Tdiv"
        png_out = os.path.join(output_folder, "Num shed genes vs over time.png")
        fig = plt.figure(figsize=(4, 4), dpi=100)
        label = "Num shed genes vs over time"
        specKS_shed_genes = [0,3333-3326,3333-3296,3333-3260]
        demographiKS_shed_genes = [0, 3333-3326, 3333-3296, 3333-3260]
        t_wgd = [10000, 100000, 500000, 1000000]
        avg_WGD_gene_lifespan = 31e6 / math.log(2) #half life to mean life conversion

        t_test = [10000, 50000, 100000, 500000, 1000000]
        shed_fraction = [1.0 - avg_WGD_gene_lifespan*expon.pdf(
            t, loc=0, scale=avg_WGD_gene_lifespan) for t in t_test]
        expectations=[f*3333 for f in shed_fraction]

        for i in range(0, len(specKS_shed_genes)):
            plt.scatter(t_wgd[i],specKS_shed_genes[i], c='r', marker="x", label="SpecKS")
            plt.scatter(t_wgd[i],demographiKS_shed_genes[i], c='b', marker="+", label="DemographiKS")

        plt.plot(t_test, expectations, c='gray',  label="expectations")

        plt.title(label)
        plt.legend()
        # plt.ylim([0,1200])
        plt.xlabel("Time (gen)")
        plt.ylabel("Num shed genes")
        plt.savefig(png_out)
        plt.clf()
        plt.close()

        self.assertEqual(True, True)  # add assertion here

if __name__ == '__main__':
    unittest.main()


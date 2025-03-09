import unittest

from figure_generation.ks_plot_aggregations_auto_vs_allo import make_Tc_Ks_Allo_vs_Auto_fig_with_subplots


class TestAlloVsAuto_Tdiv(unittest.TestCase):

    def test_allo_vs_auto_Ks_At_like_for_varying_varying_Tdiv_times(self):


        demographiKS_auto_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto/Auto_vs_Tdiv'
        demographiKS_allo_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/Auto/Auto_vs_Tdiv'
        #demographiKS_allo_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Tdiv'


        #specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Tdiv'
        demographics_allo_run_list=[False,'Allo_Twgd100v3_m03d05y2025_h13m17s02',
                    'Allo_Twgd1000v4_m03d07y2025_h08m05s29',
                    'Allo_Twgd5000v4_m03d09y2025_h13m42s12']
        #'Allo_Twgd5000v4_m03d09y2025_h13m42s12'
        demographics_auto_run_list=[False,
                                    'Auto_Twgd100v2_m03d04y2025_h17m47s07',
                                    'Auto_Twgd1000v2_m03d05y2025_h10m18s18',
                    'Auto_Twgd5000v2_m03d05y2025_h10m18s07']



        xmax_Ks_allo = [0.1 for f in demographics_auto_run_list] #0.001  # max(demographiKS_ks_results)
        xmax_Ks_auto = [0.1 for f in demographics_auto_run_list] #0.001  # max(demographiKS_ks_results)
        xmax_Tc = [80000 for f in demographics_auto_run_list]
        #xmax_Ks = [0.0005 for f in demographics_allo_run_list]
        bin_sizes_Ks_allo = [xmax_Ks_i/50.0 for xmax_Ks_i in xmax_Ks_allo]
        bin_sizes_Ks_auto = [xmax_Ks_i / 50.0 for xmax_Ks_i in xmax_Ks_auto]
        bin_sizes_Tc = [xmax_Tc_i / 50.0 for xmax_Tc_i in xmax_Tc ]

        run_list_name="Ks_for_Allo_and_Auto_At_like_varying_varying_Tdiv"
        ymax_KS = [150 for f in demographics_auto_run_list]
        ymax_Tc = [False for f in demographics_allo_run_list]
        #show_KS_predictions=[True,True,True]
        show_KS_predictions = [False, False, False]
        suptitle = "Allo and Auto Ks histograms\n" + \
                                  "Recombination rate = 8e-9, Ne and BI constant"


        bin_sizes_Ks_array=[bin_sizes_Ks_auto,bin_sizes_Ks_allo]
        xmax_Ks_array=[xmax_Ks_auto,xmax_Ks_allo]
        ymax_Ks_array=[ymax_KS,ymax_KS]

        make_Tc_Ks_Allo_vs_Auto_fig_with_subplots(bin_sizes_Ks_array, bin_sizes_Tc,
                                          demographiKS_allo_out_path, demographics_allo_run_list, run_list_name,
                                          demographics_auto_run_list, demographiKS_auto_out_path,
                                     xmax_Ks_array, xmax_Tc, ymax_Ks_array,  ymax_Tc,suptitle,
                                     show_KS_predictions)

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

        bin_sizes_Tc = [200,200, 200, 200,200]

        xmax_Ks_allo = [0.025,0.025,0.025,0.025,0.025] #0.001  # max(demographiKS_ks_results)
        xmax_Ks_auto = [0.01,0.01,0.01,0.01,0.01] #0.001  # max(demographiKS_ks_results)
        xmax_Tc = [False,False,False,False,False]
        #xmax_Ks = [0.0005 for f in demographics_allo_run_list]
        bin_sizes_Ks_allo = [xmax_Ks_i/100.0 for xmax_Ks_i in xmax_Ks_allo]
        bin_sizes_Ks_auto = [xmax_Ks_i / 100.0 for xmax_Ks_i in xmax_Ks_auto]

        run_list_name="Ks_for_Allo_and_Auto_varying_varying_Tdiv"
        ymax_KS = [800,800,800,800,800]
        ymax_Tc = [False for f in demographics_allo_run_list]
        #show_KS_predictions=[True,True,True]
        show_KS_predictions = [False, False, False]
        suptitle = "Allo and Auto Ks histograms\n" + \
                                  "Recombination rate = 8e-9, Ne and BI constant"


        bin_sizes_Ks_array=[bin_sizes_Ks_auto,bin_sizes_Ks_allo]
        xmax_Ks_array=[xmax_Ks_auto,xmax_Ks_allo]
        ymax_Ks_array=[ymax_KS,ymax_KS]

        make_Tc_Ks_Allo_vs_Auto_fig_with_subplots(bin_sizes_Ks_array, bin_sizes_Tc,
                                          demographiKS_allo_out_path, demographics_allo_run_list, run_list_name,
                                          demographics_auto_run_list, demographiKS_auto_out_path,
                                     xmax_Ks_array, xmax_Tc, ymax_Ks_array,  ymax_Tc,suptitle,
                                     show_KS_predictions)

        self.assertEqual(True, True)  # add assertion here

if __name__ == '__main__':
    unittest.main()


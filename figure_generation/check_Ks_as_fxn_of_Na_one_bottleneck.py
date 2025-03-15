import unittest

from figure_generation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots


class TestKsByNa(unittest.TestCase):

    def test_Ks_for_varying_Na_one_bottleneck_Nb_fixed_at_100(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_Na/Nb100'
        specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'

        demographics_run_list = [False,
                                 'KSvsNa_100Na_100Nb_m03d12y2025_h10m32s13',
                                 'KSvsNa_500Na_100Nb_m03d12y2025_h10m32s17',
                                   'KSvsNa_1000Na_100Nb_m03d12y2025_h10m30s08',
                                 'KSvsNa_5000Na_100Nb_m03d12y2025_h10m29s20'
                                  ]
        specks_TE5_run_list = [False,False,False,False,False,False,False]


        xmax_Ks = [0.1,0.1,0.1,0.1,0.4 ]
        bin_sizes_Ks = [xmax_KS_i/25 for xmax_KS_i in xmax_Ks]

        #xmax_Tc = [5000 for f in demographics_run_list ]
        xmax_Tc = [2000,2000,4000,8000,40000]
        bin_sizes_Tc = [xmax_Tc_i/25 for xmax_Tc_i in xmax_Tc]


        ymax_KS = [False for f in demographics_run_list]
        ymax_Tc = [False for f in demographics_run_list]

        run_list_name = "Ks_for_varying_Na_constantNb100"
        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate

        show_KS_predictions = [False, False, False]
        suptitle = "DemographiKS Ks histograms\n"
        include_annotation=False
        plot_title_lamda = lambda config: "Ks at Tnow\n"+ "Na:" + str(config.ancestral_Ne)
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, specks_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions, include_annotation, plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here


    def test_Ks_for_varying_Na_one_bottleneck_Nb_fixed_at_10K(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/KS_vs_Na/Nb10K'
        specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'

        #full, w/Ne 10K
        demographics_run_list = [False,
                                  'KSvsNa_100Na_10KNb_m03d13y2025_h09m46s52',
                                 'KSvsNa_500Na_10KNb_m03d13y2025_h09m46s52',
                                    'KSvsNa_1000Na_10KNb_m03d12y2025_h10m38s39',
                                 'KSvsNa_5000Na_10KNb_m03d12y2025_h10m38s41']
        specks_TE5_run_list = [False,False,False,False,False,False,False]


        xmax_Ks = [0.8 for f in demographics_run_list ]
        bin_sizes_Ks = [xmax_KS_i/50 for xmax_KS_i in xmax_Ks]

        xmax_Tc = [80000 for f in demographics_run_list ]
        bin_sizes_Tc = [xmax_Tc_i/50 for xmax_Tc_i in xmax_Tc]


        ymax_KS = [False for f in demographics_run_list]
        ymax_Tc = [False for f in demographics_run_list]

        run_list_name = "Ks_for_varying_Na_constant_Nb10K_figR-NaNb4"
        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate

        show_KS_predictions = [False, False, False]
        suptitle = "SLiM and SpecKS Ks histograms\n"
        include_annotation=False
        plot_title_lamda = lambda config: "Ks at Tnow\n"+ "Na:" + str(config.ancestral_Ne)
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_run_list, run_list_name,
                                     specks_TE5_run_list, specks_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, ymax_Tc,
                                     suptitle, show_KS_predictions,include_annotation, plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()

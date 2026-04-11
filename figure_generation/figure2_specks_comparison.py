import unittest
from figure_generation.multi_row_Ks_aggregation_plots import make_multi_row_Ks_fig_with_subplots


class MyTestCase(unittest.TestCase):

    def test_fig2_DGKS_vs_SPECKSKs(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Ne_v2'
        specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Ne_v2'

        r1_demographics_TE_run_list = [
                                    'DGKS_Ne10_Fig1row1_v2_m03d16y2026_h17m56s18',
                                    'DGKS_Ne50_Fig1row1_v2_m03d27y2026_h10m28s27',
                                    'DGKS_Ne100_Fig1row1_v2_m03d16y2026_h17m56s20',
                                    'DGKS_Ne500_Fig1row1_v2_m03d18y2026_h10m11s04',
                                    'DGKS_Ne1000_Fig1row1_v5_m04d01y2026_h14m06s03']


        r1_specks_TE_run_list = [
        'specks_TE10_m03d09y2026_h12m56s47','specks_TE50_m03d27y2026_h12m41s14',
                              'specks_TE100_m03d09y2026_h13m28s52',
                              'specks_TE500_m03d09y2026_h14m00s55',
                              'specks_TE1000_m03d09y2026_h13m30s53',
                              'specks_TE5000_m03d09y2026_h14m00s57']


        xmax_Ks = [0.05 for f in r1_demographics_TE_run_list]
        #xmax_Ks = [0.10 for f in demographics_TE_run_list]

        bin_sizes_Ks = [xmax_Ks_i / 50 for xmax_Ks_i in xmax_Ks]
        xmax_Tc = [1000,1000,2000,2000,10000,20000,40000, 80000]
        bin_sizes_Tc =[xmax_Tc_i / 50 for xmax_Tc_i in xmax_Tc]
        ymax_Tc = [False for f in r1_demographics_TE_run_list]
        run_list_num = "DGKS_1000_gen_by_Ne_fast_mut_rate_Fig R-Ne1.multirow."
        ymax_Ks = [200 for f in r1_demographics_TE_run_list ]
        #ymax_Ks = [600 for f in demographics_TE_run_list]

        suptitle = "DGKS vs SpecKS, Tcoal and Ks"
        show_KS_predictions=[False,False,False]
        include_annotation=False
        plot_title_lamda = lambda config: "Ks at Tnow\n"+ "Ne:" + str(config.ancestral_Ne)
        #which_plot_panels_to_show_legend = [1,2,3,4]
        which_plot_panels_to_show_legend = []
        make_multi_row_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, r1_demographics_TE_run_list, run_list_num,
                                     r1_specks_TE_run_list, specks_out_path,
                                     xmax_Ks, xmax_Tc, ymax_Ks, ymax_Tc,
                                      suptitle, show_KS_predictions,
                                     include_annotation,which_plot_panels_to_show_legend,plot_title_lamda)

        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()

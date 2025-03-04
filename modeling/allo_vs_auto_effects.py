import math
import os
import unittest
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import config
import ks_parsers
from figure_generation import curve_fitting


class MyTestCase(unittest.TestCase):

    MBE_dpi=350
    def test_coffee_histogram(self):

        coffee_num=20#coffee_num=17
        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"
        hist_comparison_out_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_vs_allo_testing"
        specks_csv_file = "Allo_Coffea{0}_ML_rep0_Ks_by_GeneTree.csv".format(coffee_num)
        ksrates_csv_file="coffea.ks.tsv"

        splat=specks_csv_file.split("_")
        species_run_name=splat[0]+splat[1]
        real_full_path=os.path.join(ksrates_out_folder,ksrates_csv_file)

        #bin_size=0.001
        bin_size=0.002
        max_Ks=0.2
        color='blue'
        wgd_ks=0.015
        density = False

        shape=0.5 #sigma
        scale = 0.02#two_Ne#going from scale 1 to 2 halves the height and makes it twice as wide
        loc = 0.01 #predic loc=0 when the mean is over the scale value.
        amp= 1.0  #0.5*num_genes*two_Ne*bin_size
        p0 = [amp,shape,loc,scale]
        lognorm_fit_range=[0.001,0.2]
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        if not os.path.exists(hist_comparison_out_folder):
            os.makedirs(hist_comparison_out_folder)

        out_png1 = os.path.join(hist_comparison_out_folder, "real_" + species_run_name + "_out.png")
        make_simple_histogram(real_ks_results, p0,lognorm_fit_range,
                              species_run_name, bin_size, color, wgd_ks,
                                           max_Ks, density, out_png1)

    def test_poplar_histogram(self):

        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"
        hist_comparison_out_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_vs_allo_testing"
        ksrates_csv_file = "poplar.ks.tsv"
        species_run_name ='Populus trichocarpa'
        real_full_path = os.path.join(ksrates_out_folder, ksrates_csv_file)
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        bin_size = 0.002
        max_Ks = 0.5
        color = 'blue'
        wgd_ks=0.21
        density = False

        shape=0.5 #sigma
        scale = 0.5#two_Ne#going from scale 1 to 2 halves the height and makes it twice as wide
        loc = 0.02 #predic loc=0 when the mean is over the scale value.
        amp= 200  #0.5*num_genes*two_Ne*bin_size
        p0 = [amp,shape,loc,scale]
        lognorm_fit_range=[0.12,max_Ks]
        out_png1 = os.path.join(hist_comparison_out_folder, "real_" + species_run_name + "_out.png")
        make_simple_histogram(real_ks_results, p0,lognorm_fit_range,
                              species_run_name, bin_size, color, wgd_ks,
                                           max_Ks, density, out_png1)

    def test_maize_histogram(self):

        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"
        hist_comparison_out_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/Auto_vs_allo_testing"
        ksrates_csv_file="mays.ks.tsv"
        species_run_name='Zea mays'
        real_full_path=os.path.join(ksrates_out_folder,ksrates_csv_file)
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        bin_size=0.002
        max_Ks=0.3
        color='blue'
        wgd_ks=0.13
        density = False

        shape=1.8 #sigma
        scale = 0.16#two_Ne#going from scale 1 to 2 halves the height and makes it twice as wide
        loc = 0.05 #predic loc=0 when the mean is over the scale value.
        amp= 300  #0.5*num_genes*two_Ne*bin_size
        lognorm_fit_range=[0.10,max_Ks]
        p0 = [amp,shape,loc,scale]

        out_png1 = os.path.join(hist_comparison_out_folder, "real_" + species_run_name + "_out.png")
        make_simple_histogram(real_ks_results, p0,lognorm_fit_range,
                              species_run_name, bin_size, color, wgd_ks,
                                           max_Ks, density, out_png1)

def make_simple_histogram(Ks_results, p0, lognorm_fit_range,
                          species_name, bin_size, color,WGD_ks, max_Ks, density, out_png):

    # MBE says: 600 - 1200 dpi for line drawings
    # and 350 dpi for color and half-tone artwork)
    fig = plt.figure(figsize=(10, 10), dpi=MyTestCase.MBE_dpi)
    x = Ks_results
    # print(PAML_hist_out_file)
    label="hist for " + os.path.basename(out_png).replace("_out.png","")
    if max_Ks:
        bins = np.arange(bin_size, max_Ks + 0.1, bin_size)
        hist_ys, bins, patches = plt.hist(x, bins=bins, facecolor=color, alpha=0.25,
                                    label=label, density=density)
        plt.xlim([0, max_Ks * (1.1)])

    bins_centers = [0.5 * (bins[i] + bins[i + 1]) for i in range(0, len(bins) - 1)]
    fit_xs=[]
    fit_ys=[]
    for i in range(0,len(bins_centers)):
        x_value=bins_centers[i]
        if x_value >= lognorm_fit_range[0]:
            if x_value < lognorm_fit_range[1]:
                fit_xs.append(x_value )
                fit_ys.append(hist_ys[i])
    fit_curve_ys_ln2, xs_for_wgd, popt = \
        curve_fitting.fit_curve_to_xs_and_ys(fit_xs,fit_ys, curve_fitting.wgd_lognorm2, p0=p0)

    if fit_curve_ys_ln2:
        popt_in_sci_notation = ["{:.2E}".format(p) for p in popt]
        plot_label = "fit popt:" + ",".join(popt_in_sci_notation)
        plt.plot(xs_for_wgd, fit_curve_ys_ln2, color='g', alpha=0.95, label=plot_label,linestyle=":")

    plt.axvline(x=WGD_ks, color='b', linestyle='-', label="WGD paralog start")
    num_pairs=sum(hist_ys)
    num_after_wgd=sum([hist_ys[i] for i in range(0,len(hist_ys)) if bins[i] > WGD_ks])
    plt.legend()
    plt.xlabel("Ks")
    plt.ylabel("Count in Bin")
    plt.title("Ks histogram for {0}.\n{1} pairs of genes. ~{2} retained from WGD.".format(
        species_name,num_pairs,num_after_wgd))
    plt.savefig(out_png)
    plt.clf()
    plt.close()

    return [hist_ys,bins]

if __name__ == '__main__':
    unittest.main()
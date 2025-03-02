import math
import os
import unittest
import random

from matplotlib.lines import lineStyles
from scipy.optimize import curve_fit
from scipy.stats import lognorm, norm, chisquare, pearsonr
import numpy as np
from matplotlib import pyplot as plt
from mkl import second
from scipy.stats import expon

from figure_generation import curve_fitting
from figure_generation.colors_for_figures import lighten_color
from figure_generation.curve_fitting import fit_curve_to_xs_and_ys, lognorm_by_sigma_mu, wgd_lognorm2


class TestModelEffectsOfCLT(unittest.TestCase):


    def test_fractional_RC(self):
        output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/CLT_testing"

        # make a bunch of exponential distributions
        num_genes = 50000
        N = 1000
        xmax = 10000
        num_bins = 200
        bin_size = xmax / num_bins
        bins = np.arange(0, xmax, bin_size)
        bins_centers = [0.5 * (bins[i] + bins[i + 1]) for i in range(0, len(bins) - 1)]
        two_Ne = 2.0 * N
        print("Tc plot bin_size_in_time: " + str(bin_size))
        kingman = [min(num_genes,
                       (bin_size * num_genes / two_Ne) * math.e ** ((-1 * i) / two_Ne))
                   for i in bins_centers]

        CM=2*N
        histograms_from_resampled_distributions_ys = [kingman]
        data_labels = ["kingman"]#, "some_RC_loss"]
        bin_size=bins_centers[1]-bins_centers[0]
        for percent_RC in [0.1]:

            some_RC_loss = [kingman[b]*(1-percent_RC) for b in range(0,len(bins_centers))]
            some_RC_gain = []
            for b in range(0,len(bins_centers)):
                x=bins_centers[b]
                if x < N:
                    some_RC_gain.append(0)
                else:
                    x0=2*x - CM
                    p0= (bin_size * num_genes / two_Ne) * math.e ** ((-1 * x0) / two_Ne)
                    some_RC_gain.append(percent_RC*p0)

            #looks like exp
            #gain2 = [bin_size*wgd_lognorm2(i,num_genes,1,-CM,CM*2)
            #           for i in bins_centers]

            #close but too far right
            #wgd_lognorm2(x, amp, shape, loc, scale)
            #gain2 = [bin_size*wgd_lognorm2(i,num_genes,1,CM/2,CM*2)
            #           for i in bins_centers]

            shape=1/4
            loc=CM/2-2000-1500
            gain2 = [bin_size*wgd_lognorm2(i,num_genes,shape,loc,CM*2)
                       for i in bins_centers]

            #new_dist=[some_RC_loss[b]+gain2[b]  for b in range(0,len(bins_centers))]
            #max_value=max(new_dist)
            #max_bin_index = new_dist.index(max_value)
            #maxx = bins_centers[max_bin_index]

            #data_labels = data_labels + ["gain_RC"+str(percent_RC),"new_dist, xmax at " + str(maxx)]
            #histograms_from_resampled_distributions_ys = histograms_from_resampled_distributions_ys + \
            #        [gain2, new_dist]

            param_string = str(bin_size) + " " + str(num_genes) + " " + str(shape) + " " + str(loc) + " " + str(CM * 2)

            data_labels = data_labels + ["gain_RC"+str(percent_RC) + "PARAMS: " + param_string
                                         ]
            histograms_from_resampled_distributions_ys = histograms_from_resampled_distributions_ys + \
                    [gain2]

        png_out_i = os.path.join(output_folder, "af_Kingman_CLT_0RC_" + param_string + ".jpg")
        plot_mrca_x_y(histograms_from_resampled_distributions_ys,
                      bins_centers, data_labels, png_out_i)
        csv_out = os.path.join(output_folder, "af_hist_data.csv")
        write_histogram_distribution_data(csv_out, data_labels, bins_centers, histograms_from_resampled_distributions_ys)

    # "A function that can be generalized to both a Gaussian (normal distribution) "
    # "and an exponential distribution is the ")generalized normal distribution"
    # (also known as the "exponential power distribution"); by adjusting a single parameter
    # within this function, you can smoothly transition between the characteristics
    # of a Gaussian and an exponential distribution depending on the parameter value"
    # https://en.wikipedia.org/wiki/Generalized_normal_distribution
    # https://en.wikipedia.org/wiki/Laplace_transform
    def test_effects_of_CLT(self):

        output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/CLT_testing"

        #make a bunch of exponential distributions
        num_genes=50000
        N=1000
        xmax=10000
        num_bins=200
        bin_size = xmax/num_bins
        bins = np.arange(0, xmax, bin_size)
        bins_centers = [0.5 * (bins[i] + bins[i + 1]) for i in range(0, len(bins) - 1)]
        list_of_expon_distributions=[]
        for num_combined_distributions in range(0,100):
            random.seed(num_combined_distributions*10+10)
            data=list(theoretical_coalescent(num_genes, N))
            list_of_expon_distributions.append(data)


        #build a new Tc drawing from multiple different distributions
        histograms_from_resampled_distributions_ys=[]
        fit_distributions_ys=[]
        data_labels=[]


        two_Ne = 2.0 * N
        print("Tc plot bin_size_in_time: " + str(bin_size))
        kingman = [min(num_genes,
                       (bin_size * num_genes / two_Ne) * math.e ** ((-1 * i) / two_Ne))
                   for i in bins[0:-1]]
        histograms_from_resampled_distributions_ys.append(kingman)
        data_labels.append("Perfect Kingman distribution")


        bin_size = bins[1] - bins[0]
        #wgd_normal(x, amp, mu, sig):
        #popt = [num_genes * bin_size, two_Ne- (bin_size/2),bin_size]
        #by CTL, sample sigma = true sigma / sqrt(#samples)
        #(picked "80" because thats the biggest #RC we tested)
        sample_sigma = two_Ne / math.sqrt(80)
        popt = [num_genes * bin_size, two_Ne, sample_sigma]
        #gaussian_ys=[curve_fitting.wgd_normal(x, *popt) for x in bins[0:-1]]
        gaussian_ys = [curve_fitting.wgd_normal(x, *popt) for x in bins_centers]
        histograms_from_resampled_distributions_ys.append(gaussian_ys)
        data_labels.append("Perfect Gaussian distribution")

        new_data = [d for d in list_of_expon_distributions[0]]
        png_out_i = os.path.join(output_folder, "CLT_0RC.png")
        dgx_hist_ys, bins = plot_mrca_histogram(new_data, bins, png_out_i)
        histograms_from_resampled_distributions_ys.append(dgx_hist_ys)
        data_labels.append("0 RC events simulated")

        data_colors = ['darkblue','darkred', 'black']

        #list_of_fraction_RC_events_to_test = [0.05,0.1,0.2,0.3,0.5,0.75]
        #blue_lightening=[0.75,0.70,0.65,0.6,0.4,0.3]
        list_of_fraction_RC_events_to_test = []
        blue_lightening=[]
        for fraction in list_of_fraction_RC_events_to_test:

            new_data = [d for d in list_of_expon_distributions[0]]
            for k in range(0, int(num_genes*fraction)):
                Tc1=list_of_expon_distributions[0][k]
                Tc2=list_of_expon_distributions[1][k]
                new_data[k] = 0.5 * (Tc1 + Tc2)

            png_out_i = os.path.join(output_folder, "CLT_" + str(fraction) + "%RC.png")
            dgx_hist_ys, bins = plot_mrca_histogram(new_data, bins, png_out_i)
            histograms_from_resampled_distributions_ys.append(dgx_hist_ys)
            data_labels.append(str(fraction) + " %RC events simulated")

        new_colors=[lighten_color('blue', amount=a) for a in blue_lightening]
        data_colors= data_colors+ new_colors

        #list_of_num_RC_events_to_test=[1,2,5,10,20,50,80]
        #red_lightening=[0.2,0.5,0.6,0.65,0.7,0.75]
        list_of_num_RC_events_to_test=[1]
        red_lightening=[]
        for num_RC_events in list_of_num_RC_events_to_test:
            num_combined_distributions=num_RC_events+1
            new_data=[0 for k in range(0, num_genes)]
            for k in range(0, num_genes):
                 Tcs_for_Kth_genes=[list_of_expon_distributions[j][k]
                                    for j in range(0, num_combined_distributions)]
                 avg_Tc_for_Kth_gene= sum(Tcs_for_Kth_genes)/num_combined_distributions
                 new_data[k]=avg_Tc_for_Kth_gene#
            png_out_i = os.path.join(output_folder, "CLT_" + str(num_RC_events) +"RC.png")
            dgx_hist_ys, bins_xs = plot_mrca_histogram(new_data, bins, png_out_i)
            histograms_from_resampled_distributions_ys.append(dgx_hist_ys)

            data_labels.append(str(num_RC_events) + " RC events simulated")

        data_colors.append("gray")
        new_colors=[lighten_color('red', amount=a) for a in red_lightening]
        data_colors= data_colors+ new_colors
        #plot the resampled distributions
        csv_out = os.path.join(output_folder, "hist_data_from_resampled_CLT_RC.csv")
        write_histogram_distribution_data(csv_out, data_labels, bins, histograms_from_resampled_distributions_ys)

        shape=1  #sigma
        scale = two_Ne#going from scale 1 to 2 halves the height and makes it twice as wide
        loc = 0 #predic loc=0 when the mean is over the scale value.
        amp=0.5*num_genes*two_Ne*bin_size
        p2 = [amp,shape,loc,scale]


        png_out = os.path.join(output_folder, "composite_CLT_RC_2.png")
        csv_out = os.path.join(output_folder, "fits_CLT_RC.csv")
        histograms_from_fit_distributions_ys = plot_composite_tc_2(histograms_from_resampled_distributions_ys,p2,
                          bins,data_labels,data_colors,png_out, csv_out)

        csv_out = os.path.join(output_folder, "hist_data_from_fit_CLT_RC.csv")
        write_histogram_distribution_data(csv_out, data_labels, bins, histograms_from_fit_distributions_ys)

        self.assertEqual(True, True)  # add assertion here

def plot_composite_tc_2(list_of_simulated_data, p0, bins,
                      data_labels, data_colors, png_out, csv_out):

    fig, ax = plt.subplots(1, 2, figsize=(20, 10))
    label = "Coalescent Times For Genes From Sampled From Polyploid Ohnologs"
    bins_centers = [ 0.5*(bins[i]+bins[i+1]) for i in range(0,len(bins)-1)]
    fit_distributions_popts =[]
    fit_distributions_ys =[]
    RCs=[]; true_max=[]; max_exp=[]

    for i in range(0, len(list_of_simulated_data)):

        data_label=data_labels[i]
        simulated_data= list_of_simulated_data[i]
        color=data_colors[i]
        ax[0].plot(bins_centers, simulated_data, color=color, alpha=0.95, label=data_label)

        # expected_peak
        if i > 3:
            num_RC = float(data_label.split(" ")[0])
            expected_cm=p0[-1]

            fit_curve_ys_ln2, xs_for_wgd, popt = \
                curve_fitting.fit_curve_to_xs_and_ys(bins_centers,simulated_data, curve_fitting.wgd_lognorm2, p0=p0)

            popt_in_sci_notation=["{:.2E}".format(p) for p in popt]
            plot_label= "fit popt:"  + ",".join(popt_in_sci_notation)
            ax[0].plot(xs_for_wgd, fit_curve_ys_ln2, color="k", alpha=0.95, label=plot_label,
                     linestyle=":")
            fit_distributions_popts.append(popt)
            fit_distributions_ys.append(fit_curve_ys_ln2)

            #expected vs true peak
            #find peak
            max_value=max(fit_curve_ys_ln2)
            max_index = fit_curve_ys_ln2.index(max_value)
            next_max_value=fit_curve_ys_ln2[max_index+1]
            max_x = bins_centers[max_index]

            #note, random effects will be bigger closer to the origin
            #I should set a random seed...
            if "%" in data_label:
                max_x_exp = expected_cm * num_RC / (num_RC + 1)
                max_x_exp = max(bins_centers[0],max_x_exp)
                diff_between_bin0_and_bin1=max_value-next_max_value
                num_values_affected=num_RC*50000 #from num genes
                if num_values_affected < diff_between_bin0_and_bin1:
                    max_x_exp = bins_centers[0]
            else:
                max_x_exp= expected_cm * num_RC / (num_RC + 1)
            print(str(num_RC) + " true_vs_exp:\t" + str(max_x) + "\t" + str(max_x_exp))
            RCs.append(num_RC)
            true_max.append(max_x)
            max_exp .append(max_x_exp)

            ax[1].scatter(max_x,max_x_exp, color=color, alpha=0.95)
            ax[0].scatter(max_x_exp,0, color=color, alpha=0.95)

    bin_size=bins_centers[1]-bins_centers[0]
    num_genes = 50000
    N = 1000
    CM=2*N


    # basically perfect3 <great fit with nice numbers!!
    shape = 0.5
    loc = -1000 #should be -N
    scale = 2500 #should be about 2Ne (from Kingman std dev) but booted over for lognorm
    #scale from (1000*sqrt(2*pi)). see https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.lognorm.html

    gain2 = [bin_size*wgd_lognorm2(i,num_genes,shape,loc,scale ) for i in bins_centers]

    param_string = str(bin_size) + " " + str(num_genes) + " " + str(shape) + \
                    " " + str(loc) + " " + str(scale)

    ax[0].plot(bins_centers, gain2 , color="orange", alpha=0.95, label="RC0\n"+\
         "PARAMS: " +param_string)
      #  +str(bin_size) + " " + str(num_genes) + " " + str(1) + " " + str(0) + " " + str(CM * 2))

    expt_csv = csv_out.replace(".csv", "_expectations.csv")
    write_some_data(expt_csv, RCs, true_max, max_exp)
    ax[1].plot([0,2000], [0,2000], color='gray', alpha=0.95, linestyle="-")

    plt.title(label)
    ax[0].legend()
    ax[0].set(ylim=[0,1400])
    ax[0].set(xlabel="MRCA time")
    ax[0].set(ylabel="# genes in bin")
    ax[1].set(xlabel="true position of Ks max")
    ax[1].set(ylabel="expected position of Ks max")
    ax[0].set(title=label)
    plt.savefig(png_out.replace(".png",param_string + ".png"))
    plt.clf()
    plt.close()
    return fit_distributions_ys


def write_histogram_distribution_data(csv_out, data_labels, bins, resampled_distributions_ys):
    with open(csv_out, 'w') as f:

        bins_str = [str(b) for b in bins]
        f.write("Bins\t"+ "\t".join(bins_str) +"\n")

        f.write("RC\tDistributionData\n")
        for i in range(0, len(resampled_distributions_ys)):
            data_as_str = [str(d) for d in resampled_distributions_ys[i]]
            f.write(data_labels[i] + "\t" + "\t".join(data_as_str) + "\n")

def write_some_data(csv_out, RCs, true_max,max_exp):

    with open(csv_out, 'w') as f:
        f.write("RCs\ttrue_max\tmax_exp\n")
        for i in range(0, len(RCs)):
            data=[RCs[i], true_max[i],max_exp[i]]
            data_as_str = [str(d) for d in data]
            f.write("\t".join(data_as_str)+ "\n")

def plot_mrca_histogram(theoretical_mrcas_by_gene, bins, png_out):

    fig = plt.figure(figsize=(10, 10), dpi=350)
    label = "Coalescent Times For Genes From Sampled Ohnologs"
    num_genes=len(theoretical_mrcas_by_gene)
    dgx_hist_ys, bins, patches = plt.hist(theoretical_mrcas_by_gene, bins=bins, facecolor='c', alpha=0.25,
                 label='Theoretical Tcoal by gene (total: ' + str(num_genes) + ')',
                 density=False)

    plt.title(label)
    plt.legend()
    plt.xlabel("MRCA time")
    plt.ylabel("# genes in bin")
    plt.savefig(png_out)
    plt.clf()
    plt.close()
    return dgx_hist_ys, bins

def plot_mrca_x_y(data_y_list, bins_x, data_labels, png_out):

    fig = plt.figure(figsize=(10, 10), dpi=350)
    label = "Coalescent Times For Genes From Sampled Ohnologs"

    #, color = 'c'
    for i in range(0,len(data_y_list)):
        plt.plot(bins_x, data_y_list[i], alpha=1,
                 label= data_labels[i])

    plt.title(label)
    plt.legend()
    plt.xlabel("MRCA time")
    plt.ylabel("# genes in bin")
    plt.savefig(png_out)
    plt.clf()
    plt.close()
    return

def theoretical_coalescent(num_genes,N):
    #Co.T=(1/2N)*e^-((t-1)/2N))
    loc=0
    xscale=2.0*N
    random_draws_from_distribution = expon.rvs(loc=loc, size=num_genes, scale=xscale)
    return random_draws_from_distribution

if __name__ == '__main__':
    unittest.main()

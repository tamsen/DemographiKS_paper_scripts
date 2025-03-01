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

        list_of_fraction_RC_events_to_test = [0.05,0.1,0.2,0.3,0.5,0.75]
        blue_lightening=[0.75,0.70,0.65,0.6,0.4,0.3]
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

        list_of_num_RC_events_to_test=[1,2,5,10,20,50,80]
        red_lightening=[0.2,0.5,0.6,0.65,0.7,0.75]
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

def write_histogram_distribution_data(csv_out, data_labels, bins, resampled_distributions_ys):
    with open(csv_out, 'w') as f:

        bins_str = [str(b) for b in bins]
        f.write("Bins\t"+ "\t".join(bins_str) +"\n")

        f.write("RC\tDistributionData\n")
        for i in range(0, len(resampled_distributions_ys)):
            data_as_str = [str(d) for d in resampled_distributions_ys[i]]
            f.write(data_labels[i] + "\t" + "\t".join(data_as_str) + "\n")

def write_fit_parameters(csv_out, data_labels, bin_centers, fit_distributions_popts, fit_distributions_ys):

    data_list=[]
    with open(csv_out, 'w') as f:
        bins_str = [str(b) for b in bin_centers]
        f.write("Bins\t"+ "\t".join(bins_str) +"\n")
        num_xs=len(bin_centers)
        fraction_num_xs= 1.0 /float(num_xs)
        #popt_labels=["popt"+str(i) for i in range(0,len(fit_distributions_popts[0]))]
        popt_labels=["amp", "shape", "loc", "scale",
                     "cm_x", "peak_x", "std_dev",
                     "cm_x_exp", "peak_x_exp", "std_dev_exp",
                     "sigma_exp","mu_exp",
                     "amp_exp","shape_exp","loc_exp","scale_exp"]
        #add mean,mode,variance,skew
        # predictions based on RC

        expected_cm=2000 #2Ne

        f.write("RC\t" + "\t".join(popt_labels)  + "\n")
        for i in range(0, len(fit_distributions_ys)):

            label=data_labels[i+3]
            num_RC=float(label.split(" ")[0])
            fit_distribution=fit_distributions_ys[i]
            popt=fit_distributions_popts[i]
            popt_as_str = [str(d) for d in popt]

            sum_ys=sum(fit_distribution)
            normalizer = 1.0 / sum_ys

            #find cm
            cm_terms= [fit_distribution[j]*bin_centers[j] for j in range(0,num_xs)]
            cm_x = normalizer * sum(cm_terms)

            #find peak
            max_index = fit_distribution.index(max(fit_distribution))
            max_x = bin_centers[max_index]

            #expected_peak
            #max_x_exp = expected_cm * num_RC / (num_RC + 1)
            #max_x_exp = max(bins_centers[0], max_x_exp)

            max_x_exp= expected_cm * num_RC / (num_RC + 1)

            #std dev
            #sum (x-mu)(x-mu)*P(x)
            terms=[ fit_distribution[i]*(bin_centers[i]-max_x)**2
                    for i in range(0,len(bin_centers))]
            #normalizer=1.0/(sum_ys*(num_xs-1))
            sigma_squared=normalizer* sum(terms)
            std_dev= math.sqrt(sigma_squared)

            #expected std deviation
            std_dev_exp = expected_cm / (num_RC+1)

            #https://medium.com/towards-data-science/log-normal-distribution-a-simple-explanation-7605864fb67c#:~:text=Calculate%20median%2C%20mean%2C%20mode%20&%20variance&text=How%20do%20we%20arrive%20at,the%20mean%20(see%20here).
            A = max_x_exp
            B = expected_cm
            lnA = math.log(A)
            lnB = math.log(B)
            sigma=math.sqrt(lnB-lnA)
            mu = 2.0*lnB + lnA

            amp = 10000000.0
            shape = -1.0*(max_x_exp - expected_cm) / (2.0* max_x_exp)
            #shape = (expected_cm -max_x_exp) / expected_cm
            #shape = (max_x_exp - expected_cm) / max_x_exp

            loc = -1*max_x_exp / 2.0
            loc= ((500)*(max_x_exp - expected_cm) / (expected_cm))-1000

            scale = 1000.0 + max_x_exp

            data=[cm_x,max_x,std_dev,
                  expected_cm,max_x_exp,std_dev_exp,
                  sigma,mu,amp,shape,loc,scale]
            data_as_str = [str(d) for d in data]
            print(str(num_RC) + " true_vs_exp:\t" + str(max_x) + "\t" + str(max_x_exp))
            f.write(label + "\t" + "\t".join(popt_as_str) +
                    "\t" + "\t".join(data_as_str) +
                    "\n")

            data_list.append(data)
    return data_list

def read_distribution_tests(csv_in):
    results={}
    with open(csv_in, 'r') as f:
        lines=f.readlines()
        bins_splat=lines[0].split('\t')
        bins_as_floats = [float(b) for b in bins_splat[1:-1]]
        for l in lines[2:-1]:
            splat=l.split('\t')
            label=splat[0]
            data=splat[1:-1]
            data_as_float=[float(d) for d in data]
            results[label]=data_as_float
    return results,bins_as_floats

def plot_composite_tc_2(list_of_simulated_data, p0, bins,
                      data_labels, data_colors, png_out, csv_out):

    fig, ax = plt.subplots(1, 2, figsize=(20, 10))
    label = "Coalescent Times For Genes From Sampled From Polyploid Ohnologs"
    bins_centers = [ 0.5*(bins[i]+bins[i+1]) for i in range(0,len(bins)-1)]
    fit_distributions_popts =[]
    fit_distributions_ys =[]

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
                #max_x_exp = expected_cm * num_RC * ( num_RC / (1+num_RC))**(1/2)
                #max_x_exp = expected_cm *   num_RC**2 / (num_RC**2 + 1)
                #max_x_exp = expected_cm * num_RC / (num_RC + 1)
            else:
                max_x_exp= expected_cm * num_RC / (num_RC + 1)
            print(str(num_RC) + " true_vs_exp:\t" + str(max_x) + "\t" + str(max_x_exp))
            ax[1].scatter(max_x,max_x_exp, color=color, alpha=0.95)
            ax[0].scatter(max_x_exp,0, color=color, alpha=0.95)

    ax[1].plot([0,2000], [0,2000], color='gray', alpha=0.95, linestyle="-")
    #write_histogram_distribution_data(csv_out, data_labels, bins, fit_distributions_ys)
    #metrics_csv=csv_out.replace(".csv","_metrics.csv")
    #data_list = write_fit_parameters(metrics_csv, data_labels, bins_centers,
    #                                 fit_distributions_popts, fit_distributions_ys)


    plt.title(label)
    ax[0].legend()
    ax[0].set(xlabel="MRCA time")
    ax[0].set(ylabel="# genes in bin")
    ax[1].set(xlabel="true position of Ks max")
    ax[1].set(ylabel="expected position of Ks max")
    ax[0].set(title=label)
    plt.savefig(png_out)
    plt.clf()
    plt.close()
    return fit_distributions_ys

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

def theoretical_coalescent(num_genes,N):
    #Co.T=(1/2N)*e^-((t-1)/2N))
    loc=0
    xscale=2.0*N
    random_draws_from_distribution = expon.rvs(loc=loc, size=num_genes, scale=xscale)
    return random_draws_from_distribution

if __name__ == '__main__':
    unittest.main()

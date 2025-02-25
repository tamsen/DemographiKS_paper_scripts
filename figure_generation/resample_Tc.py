import os.path
import unittest

import config
from figure_generation.coalescent_plot_aggregation import make_Tc_fig_with_subplots
from modules import trees_file_processor


class TestResampleTc(unittest.TestCase):

    def test_ResampleTc(self):

        output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/trees_file_testing"
        trees_file="allotetraploid_bottleneck_trees_at_div_32.txt"
        input_xml_file="Inb32v1.used.xml"
        full_path_to_trees=os.path.join(output_folder,trees_file)
        full_path_to_config=os.path.join(output_folder,input_xml_file)
        config_used = config.DemographiKS_config(full_path_to_config)


        # select random ancestral genomes to calculate the T coalescent:
        pairs=[[1,2],[1,5],[1,100],[1,500]]
        for pair in pairs:

            trees_file_processor.plot_coalescent(full_path_to_trees,  pair[0],pair[1],
                                                 config_used, output_folder)
        self.assertEqual(True, True)  # add assertion here



if __name__ == '__main__':
    unittest.main()

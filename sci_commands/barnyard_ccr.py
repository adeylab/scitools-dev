#!/home/users/oconnelb/miniconda2/envs/python3.7/bin/python

import pandas as pd
import sys
import warnings
import argparse

### Ignore pandas warning about how I'm slicing the dataframe in place.  It's fine because the two slice conditions are mutually exclusive.
warnings.filterwarnings("ignore")

def main():

    ### Parse command line arguments
    parser = argparse.ArgumentParser(prog = 'calculateCCR', description = 'Calculates the median cross-contamination rate for a barnyard cells file from scitools barnyard-compare.')
    parser.add_argument('-c', nargs='?', type=str, required=True, help= ".barnyard_cells.txt file")
    parser.add_argument('-r', nargs='?', type=float, default=0.33, help="Middle range to treat as putative doublets and exclude from calculation [default =0.33, or middle third]")
    options = parser.parse_args(sys.argv[1:])

    barnyard = pd.read_table(options.c, header=None, names=["Total Reads", "Unique Human Reads", "Unique Mouse Reads", "Ratio", "Species Called"])
    mouse_in_human=barnyard.loc[barnyard["Ratio"]>(1-options.r)]
    mouse_in_human["Exact Ratio"]=1-(mouse_in_human["Unique Human Reads"]/mouse_in_human["Total Reads"])
    human_in_mouse=barnyard.loc[barnyard["Ratio"]<options.r]
    human_in_mouse["Exact Ratio"]=1- (human_in_mouse["Unique Mouse Reads"]/human_in_mouse["Total Reads"])

    sys.stdout.write("Median CCR: {:4f}".format(human_in_mouse["Exact Ratio"].median()+mouse_in_human["Exact Ratio"].median()))

if __name__=='__main__':
    main()

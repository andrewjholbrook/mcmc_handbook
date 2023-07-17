# This file generate the xml files to run BEAST analyses
# Xiang Ji
import argparse, os, re, numpy

def main():
    for steps in [4,8,16,32,64]:
        with open("mcmc_handbook/code/tuning.xml","r") as infile:
            with open("mcmc_handbook/code/timing" +str(steps)  +".xml","w") as outfile:
                for s in infile:
                    outfile.write(s.replace("nSteps=\"4\"", "nSteps=\"" + str(steps) + "\"").replace("chainLength=\"100000\"", "chainLength=\"1000\""))


if __name__ == '__main__':
    main()

# python get_locations_times.py

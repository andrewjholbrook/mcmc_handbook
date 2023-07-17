# This file generate the xml files to run BEAST analyses
# Xiang Ji
import argparse, os, re, numpy

def main():
    for steps in [4,8,16,32,64]:
        for rates in [0.2,0.3,0.4,0.5,0.6,0.7,0.8]:
            with open("mcmc_handbook/code/tuning.xml","r") as infile:
                with open("mcmc_handbook/code/tuning_" +str(steps) +"_" + str(rates) +".xml","w") as outfile:
                    for s in infile:
                        outfile.write(s.replace("nSteps=\"4\"", "nSteps=\"" + str(steps) + "\"").replace("targetAcceptanceProbability=\"0.45\"", "targetAcceptanceProbability=\"" +  str(rates) + "\"").replace("tuning.log","tuning"+str(steps)+str(rates)+".log"))


if __name__ == '__main__':
    main()

# python get_locations_times.py

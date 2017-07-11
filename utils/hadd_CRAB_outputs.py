#!/usr/bin/env python
import os, re
import multiprocessing
import logging
import commands
import math, time
import sys
from ROOT import TObject, TFile, TH1, TH1F
from array import array

#LUMI        =  35867# in pb-1

# use the following lists to include/exclude samples to be merged

blacklist = []
whitelist = []

########## DO NOT TOUCH BELOW THIS POINT ##########

import argparse

parser = argparse.ArgumentParser(description='combine the CRAB outputs into one tree')
parser.add_argument('folder', help='the crab_projects folder')
args = parser.parse_args()

if not os.path.exists(os.path.expandvars(args.folder)):
    print '--- ERROR ---'
    print '  \''+args.folder+'\' path not found'
    print '  please point to the correct path to the folder containing CRAB outputs' 
    print 
    exit()

jobs = []
#names = []

def hadd_samples(name):
    if len(whitelist)>0 and not name in whitelist: return
    if len(blacklist)>0 and name in blacklist: return
#by starting from [5:] it gets rid of "crab_" string
    os.system('hadd -f '+name[5:]+'.root '+name+'/*/*.root')
    os.system('ls '+name[5:]+'.root\n')
pass

subdirs = [x for x in os.listdir(args.folder) if os.path.isdir(os.path.join(args.folder, x))]
print subdirs

os.chdir(args.folder)
for l in subdirs:
    p = multiprocessing.Process(target=hadd_samples, args=(l,))
    jobs.append(p)
    #print p.name
    p.start()

os.system('cd ..')

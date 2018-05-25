#!/usr/bin/env python
import os
#samplenames = ['SingleMuon', 'MET']
samplenames = ['SingleMuon']
##samplename = 'SingleMuon'
##samplename = 'MET'
datasets = ['Run2018A-PromptReco-v1','Run2017A-PromptReco-v2']

#datasets = ['Run2017E-PromptReco-v1','Run2017F-PromptReco-v1']
#datasets = ['Run2017C-PromptReco-v3']
########parser#######
import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-a", "--crabaction", action="store", type="string", dest="crabaction", default="")
(options, args) = parser.parse_args()
####################
for samplename in samplenames:
    for a in datasets:
        if options.crabaction=="submit":
            os.system('echo submitting... crabConfigTrig with options... General.requestName=' + samplename + '_' + a + ' Data.inputDataset=/' + samplename + '/' + a + '/MINIAOD\n')
            os.system('crab submit -c crabConfigTrig.py General.requestName=' + samplename + '_' + a + ' Data.inputDataset=/' + samplename + '/' + a + '/MINIAOD\n')
        elif options.crabaction=="status":
            os.system('echo status -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
            os.system('crab status -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
        elif options.crabaction=="resubmit":
            os.system('echo resubmit -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
            os.system('crab resubmit -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
        elif options.crabaction=="getoutput":
            if a=='Run2017B-PromptReco-v1':
                os.system('echo getoutput -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
                os.system('crab getoutput --jobids 1-250 -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
                os.system('crab getoutput --jobids 251-543 -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
            elif a =='Run2017C-PromptReco-v2':
                os.system('echo getoutput -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
                os.system('crab getoutput --jobids 1-250 -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
                os.system('crab getoutput --jobids 251-530 -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
            elif a =='Run2017C-PromptReco-v3':
                os.system('echo getoutput -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
                os.system('crab getoutput --jobids 1-380 -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
                os.system('crab getoutput --jobids 381-751 -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
            elif a =='Run2017D-PromptReco-v1':
                os.system('echo getoutput -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
                os.system('crab getoutput --jobids 1-386 -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
                os.system('crab getoutput --jobids 387-771 -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
            elif a =='Run2017E-PromptReco-v1':
                os.system('echo getoutput -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
                os.system('crab getoutput --jobids 1-407 -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
                os.system('crab getoutput --jobids 408-815 -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
                os.system('crab getoutput --jobids 815-1221 -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
            else:
                os.system('echo getoutput -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
                os.system('crab getoutput -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
        elif options.crabaction=="kill":
            os.system('echo kill -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
            os.system('crab kill -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
        elif options.crabaction=="report":
            os.system('echo report -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
            os.system('crab report -d crab_projects_trigger_1_nov_C/crab_'+samplename+'_'+a+'\n')
        else:
            print "Invalid crab action. Please type: -a submit/status/resubmit/getoutput/kill"
            exit()
os.system('echo --------------------------\n') 

#config.General.requestName = 'SingleMuon_Run2017B-PromptReco-v1'#543
#config.Data.inputDataset = '/SingleMuon/Run2017B-PromptReco-v1/MINIAOD'
#config.General.requestName = 'SingleMuon_Run2017B-PromptReco-v2'
#config.Data.inputDataset = '/SingleMuon/Run2017B-PromptReco-v2/MINIAOD'
#config.General.requestName = 'SingleMuon_Run2017C-PromptReco-v1'
#config.Data.inputDataset = '/SingleMuon/Run2017C-PromptReco-v1/MINIAOD'#530
#config.General.requestName = 'SingleMuon_Run2017C-PromptReco-v2'
#config.Data.inputDataset = '/SingleMuon/Run2017C-PromptReco-v2/MINIAOD'
#config.General.requestName = 'SingleMuon_Run2017C-PromptReco-v3'
#config.Data.inputDataset = '/SingleMuon/Run2017C-PromptReco-v3/MINIAOD'

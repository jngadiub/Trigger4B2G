from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import os
config = config()

config.General.workArea = 'crab_projects_trigger_2018_25May2018'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/ConfFile_cfg.py'
config.JobType.inputFiles = ['data']

#config.General.requestName = 'SingleMuon_Run2017A-PromptReco-v2' #no LS in golden json
#config.Data.inputDataset = '/SingleMuon/Run2017A-PromptReco-v2/MINIAOD' #no LS in golden json
#config.General.requestName = 'SingleMuon_Run2017A-PromptReco-v3' #no LS in golden json
#config.Data.inputDataset = '/SingleMuon/Run2017A-PromptReco-v3/MINIAOD' #no LS in golden json
#config.General.requestName = 'SingleMuon_Run2017B-PromptReco-v1'
#config.Data.inputDataset = '/SingleMuon/Run2017B-PromptReco-v1/MINIAOD'
#config.General.requestName = 'SingleMuon_Run2017B-PromptReco-v2'
#config.Data.inputDataset = '/SingleMuon/Run2017B-PromptReco-v2/MINIAOD'
config.General.requestName = 'SingleMuon_Run2017C-PromptReco-v1'
config.Data.inputDataset = '/SingleMuon/Run2017C-PromptReco-v1/MINIAOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'#if Data
#config.Data.splitting = 'FileBased'#if MC

#if DATA:
config.Data.lumiMask = 'data/JSON/Cert_314472-316271_13TeV_PromptReco_Collisions18_JSON.txt'
#'Cert_294927-302654_13TeV_PromptReco_Collisions17_JSON.txt'#17.85 fbinv
#'data/JSON/Cert_294927-302343_13TeV_PromptReco_Collisions17_JSON.txt'#'data/JSON/Cert_294927-301567_13TeV_PromptReco_Collisions17_JSON.txt'#golden json
#config.Data.lumiMask = 'data/JSON/json_DCSONLY.txt' #DCS only

#older
#config.Data.lumiMask = 'data/JSON/Cert_294927-297723_13TeV_PromptReco_Collisions17_JSON.txt'
#'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'#proper way to link a json file from DQM certification
config.Data.unitsPerJob = 35

#If storing in T2
#config.Data.outLFNDirBase = '/store/user/lbenato/MET_trigger'
config.Data.publication = False
#config.Data.outputDatasetTag = 'MET_trigger'

#If willing to publish on dbs
#config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.outLFNDirBase = '/store/user/lbenato/trigger_2018_25May2018'

config.Site.storageSite = 'T2_DE_DESY'#modify with a T2 where you have writing access
config.Site.blacklist   = ['T2_FR_IPHC']

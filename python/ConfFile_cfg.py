import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import os

options = VarParsing('analysis')
options.parseArguments()

process = cms.Process("Trigger")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'/store/data/Run2017A/SingleElectron/MINIAOD/PromptReco-v2/000/296/168/00000/9ED74C00-5D4C-E711-9499-02163E01A667.root'
        #'/store/data/Run2017A/SingleMuon/MINIAOD/PromptReco-v2/000/296/168/00000/084C505D-784C-E711-8140-02163E019DA4.root'
        '/store/data/Run2017B/SingleMuon/MINIAOD/PromptReco-v2/000/298/678/00000/B48DBFCD-A466-E711-A86B-02163E01A456.root'
    )
)

process.TFileService = cms.Service( "TFileService",
    fileName = cms.string('output.root' if len(options.outputFile)==0 else options.outputFile),
    closeFileFast = cms.untracked.bool(True),
)

#-----------------------#
#     DATA FLAGS        #
#-----------------------#
isData          = ('/store/data/' in process.source.fileNames[0])
isPromptReco    = ('PromptReco' in process.source.fileNames[0])

#-----------------------#
#    VERTEX FILTER      #
#-----------------------#
import RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi
process.primaryVertexFilter = cms.EDFilter('GoodVertexFilter',
    vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
    minimumNDOF = cms.uint32(4),
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)

#-----------------------#
#     GLOBAL TAG        #
#-----------------------#
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
GT = ''
if isData:
    if isPromptReco: GT = "92X_dataRun2_Prompt_v4"
    print "data 2017, PromptReco"
else:
    GT = "90X_upgrade2017_realistic_v20"

process.GlobalTag = GlobalTag(process.GlobalTag, GT)
print 'GlobalTag loaded: ', GT

#-----------------------#
#        FILTERS        #
#-----------------------#

# JSON filter
import FWCore.PythonUtilities.LumiList as LumiList
jsonName = "Cert_294927-299420_13TeV_PromptReco_Collisions17_JSON" #golden json
process.source.lumisToProcess = LumiList.LumiList(filename = 'data/JSON/'+jsonName+'.txt').getVLuminosityBlockRange()
print "JSON file loaded: ", jsonName

process.trigger = cms.EDAnalyzer('TrigAnalyzer',
    verbose = cms.bool(True),
)


process.p = cms.Path(process.trigger)

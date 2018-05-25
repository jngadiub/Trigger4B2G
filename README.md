# Trigger4B2G
First tentative super-simplified ntuplizer to run super-easy trigger studies on 2018 early data

## Notes on trigger studies
### 25 May 2018
* 2018 data are new, hence we don't have POG blessed objects or recipes yet. We are going to read them from plain miniAOD.
* Getting muon ID from miniAOD is straightforward, jet ID a bit more difficult (but done), electron ID still work in progress (old 8XX recipes don't work anymore).
* Global Tag is the proper one for PromptReco (taken from PPD twiki).
* TrigAnalyzer.cc so far calculates MET trigger efficiencies on single muon dataset.
* Latest golden JSON files available
* MET filters being tested

## Git prerequisites
git account

git environment set

## Git instructions

In your working area, first set up the CMSSW release and initialize git:
```bash
cmsrel CMSSW_10_1_4_patch1
cd CMSSW_10_1_4_patch1/src
cmsenv
git cms-init
```

Clone the repository:

```bash
cd $CMSSW_BASE/src
mkdir TrigAnalyzer
cd TrigAnalyzer
git clone -b 2018 --single-branch https://github.com/lbenato/Trigger4B2G.git
```

## Compile and Run

Compile:
```bash
cd $CMSSW_BASE/src
scram b -j 8
```

Run:
```bash
cmsRun python/ConfFile_cfg.py
```

## CRAB
There is also a simple crabConfigTrig.py configuration file.

If you want to add together CRAB output files:
```bash
python utils/hadd_CRAB_outputs.py crab_projects_trigger
#argument: name of the crab projects folder
```

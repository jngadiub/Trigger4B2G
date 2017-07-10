# Trigger4B2G
First tentative super-simplified ntuplizer to run super-easy trigger studies on 2017 early data

## Git prerequisites
git account

git environment set

## Git instructions

In your working area, first set up the CMSSW release and initialize git:
```bash
cmsrel CMSSW_9_2_2
cd CMSSW_9_2_2/src
cmsenv
git cms-init
```

Clone the repository:

```bash
cd $CMSSW_BASE/src
git clone https://github.com/lbenato/Trigger4B2G.git
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
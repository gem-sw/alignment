# alignment
GEM alignment studies

```
scram p -n DirectoryName CMSSW_12_0_3
cd DirectoryName/src/
git clone git@github.com:gem-sw/alignment.git
cd alignment/GEM_Alignment/
cmsenv
scram b -j8
cd test
voms-proxy-init --rfc --voms cms
cmsRun analyser.py
```

Currently set to run all 3 propagation methods on a test file from 2021 Comissioning ZeroBias0 dataset

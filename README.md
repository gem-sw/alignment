# alignment
GEM alignment studies

scram p -n DirectoryName CMSSW_12_0_3
cd GEMCSCBendingAnalyzer
scram b -j8
cd test
cmsRun analyser.py

Currently set to run all 3 propagation methods on a test file from 2021 Comissioning ZeroBias0 dataset

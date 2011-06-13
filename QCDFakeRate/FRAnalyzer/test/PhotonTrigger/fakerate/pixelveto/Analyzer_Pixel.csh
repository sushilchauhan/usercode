#!/bin/tcsh
cd /uscms_data/d2/sushil/CMSSW/MonoPhoton/CMSSW_4_2_3/src
source /uscmst1/prod/sw/cms/setup/cshrc prod
cmsenv
cd /uscms_data/d2/sushil/CMSSW/MonoPhoton/CMSSW_4_2_3/src/QCDFakeRate/FRAnalyzer/test/PhotonTrigger/fakerate/pixelveto
/uscms_data/d2/sushil/CMSSW/MonoPhoton/CMSSW_4_2_3/src/QCDFakeRate/FRAnalyzer/test/PhotonTrigger/fakerate/pixelveto/myPlot.exe

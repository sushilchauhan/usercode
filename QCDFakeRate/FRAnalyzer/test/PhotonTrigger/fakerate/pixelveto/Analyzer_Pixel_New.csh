#!/bin/tcsh
cd /uscms_data/d2/sushil/CMSSW/MonoPhoton/CMSSW_3_9_7/src
source /uscmst1/prod/sw/cms/setup/cshrc prod
setenv SCRAM_ARCH slc5_ia32_gcc434
cmsenv
cd /uscms_data/d2/sushil/CMSSW/MonoPhoton/CMSSW_3_9_7/src/QCDFakeRate/Analyser/test/PhotonTrigger/fakerate/pixelveto
/uscms_data/d2/sushil/CMSSW/MonoPhoton/CMSSW_3_9_7/src/QCDFakeRate/Analyser/test/PhotonTrigger/fakerate/pixelveto/myPlot_New.exe

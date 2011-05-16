#!/bin/tcsh

setenv pwd $PWD

cat>Job_${1}.csh<<EOF
#!/bin/tcsh
#setenv X509_USER_PROXY /uscms/home/sushil/x509up_u12124 
source /uscmst1/prod/sw/cms/setup/cshrc prod
cd /uscms_data/d2/sushil/CMSSW/MonoPhoton/CMSSW_3_10_0/src
cmsenv
cd ${pwd}
#cmsDriver.py Configuration/GenProduction/python/${1}_cff.py -s GEN --condition DESIGN310_V3::All --beamspot Realistic8TeVCollision --datatier GEN-SIM --eventcontent RAWSIM --no_exec -n 10000
cmsDriver.py Configuration/GenProduction/python/${1}_cff.py -s GEN,SIM --condition DESIGN310_V3::All --beamspot Realistic8TeVCollision --datatier GEN-SIM --eventcontent RAWSIM --no_exec -n 20
cmsRun -p /uscms_data/d2/sushil/CMSSW/MonoPhoton/CMSSW_3_10_0/src/${1}_cff_py_GEN_SIM.py >  /uscms_data/d2/sushil/CMSSW/MonoPhoton/CMSSW_3_10_0/src/TimeSize_${1}
EOF

chmod 775 ${pwd}/Job_${1}.csh

cat>condor_${1}<<EOF
universe = vanilla
Executable = ${pwd}/Job_${1}.csh
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
SSC 
Error
Log 
notify_user = sushil@fnal.gov
Queue 1
EOF

sed -e 's|SSC| Output = $(Cluster)_$(Process).stdout|'   ${pwd}/condor_${1} > ${pwd}/condor_${1}_tmp
mv ${pwd}/condor_${1}_tmp ${pwd}/condor_${1}
sed -e 's|Error| Error = $(Cluster)_$(Process).stderr|'  ${pwd}/condor_${1} > ${pwd}/condor_${1}_tmp
mv ${pwd}/condor_${1}_tmp ${pwd}/condor_${1}
sed -e 's|Log| Log = $(Cluster)_$(Process).log|'   ${pwd}/condor_${1} > ${pwd}/condor_${1}_tmp
mv ${pwd}/condor_${1}_tmp ${pwd}/condor_${1}

condor_submit condor_${1}


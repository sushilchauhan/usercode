#!/bin/tcsh
setenv pwd $PWD
cp ~/tmp/x509up_u12124 /uscms/home/sushil/ 
cat>Job_${1}.csh<<EOF
#!/bin/tcsh
setenv X509_USER_PROXY /uscms/home/sushil/x509up_u12124
source /uscmst1/prod/sw/cms/setup/cshrc prod
cd /uscms_data/d2/sushil/CMSSW/MonoPhoton/CMSSW_3_9_7/src 
cmsenv
cd ${pwd}
/opt/d-cache/srm/bin/srmcp "file://localhost/${PWD}/${1}.root" "srm://cmssrm.fnal.gov:8443/11/store/user/sushil/MonoPhoton/397_Ntuples_V26/ShortRootFiles/${1}.root"
echo "Copied file"
EOF

chmod 775 $PWD/Job_${1}.csh

cat>condor_${1}<<EOF
universe = vanilla
Executable = $PWD/Job_${1}.csh
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
SSC 
Error
Log 
notify_user = schauhan@cern.ch
Queue 1
EOF

sed -e 's|SSC| Output = $(Cluster)_$(Process).stdout|'   $PWD/condor_${1} > $PWD/condor_${1}_tmp
mv $PWD/condor_${1}_tmp $PWD/condor_${1}
sed -e 's|Error| Error = $(Cluster)_$(Process).stderr|'  $PWD/condor_${1} > $PWD/condor_${1}_tmp
mv $PWD/condor_${1}_tmp $PWD/condor_${1}
sed -e 's|Log| Log = $(Cluster)_$(Process).log|'   $PWD/condor_${1} > $PWD/condor_${1}_tmp
mv $PWD/condor_${1}_tmp $PWD/condor_${1}

condor_submit condor_${1}


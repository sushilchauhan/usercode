#!/bin/tcsh
#setenv X509_USER_PROXY /uscms/home/sandhya/x509up_u45353
setenv pwd $PWD
echo ${pwd}

cat>Job_${1}.csh<<EOF
#!/bin/tcsh
source /uscmst1/prod/sw/cms/setup/cshrc prod
cd /uscms_data/d2/sushil/CMSSW/MonoPhoton/CMSSW_3_9_7/src 
cmsenv
cd ${pwd}
${pwd}/${1}.exe
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


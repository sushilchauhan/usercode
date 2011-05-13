#!/bin/tcsh


##-------GIVE INPUT HERE---------------------------------------------------
#Destination directory in dCache starting with usernaem e.g., /sushil/..
setenv destinationDir /sushil/Test 


##the directory where root files are located (directory should have only files to be transferred)
setenv sourceDir /uscmst1b_scratch/lpc1/3DayLifetime/sushil/Skim_A2

#you know this from voms-proxy-info(use yours)
setenv proxy x509up_u12124 
#-------------------------------------------------------------------------








echo '-------------------------------'
echo 'Have you done voms-proxy-init  '
echo 'Have you put your x509up       '
echo '-------------------------------'



setenv input_files inputFiles.list
setenv cmsswDir /uscms_data/d2/sushil/CMSSW_4_1_4/src

setenv pwd $PWD
cp /tmp/${proxy} /uscms/home/${USER}/ 


#create the .list file from sourceiDirec
ls -tr --format=single ${sourceDir} > inputFiles.list 


#Read the input files
set file_length=`wc -l ${input_files} | cut -c1-2`                  # Count Lines
          
 
set p = 0  
           
#run till the last line of the input files
while ($p != $file_length)
@ p = ${p} + 1
set rootFileName=`tail -n +$p ${input_files} | head -1`                  # Get name of dataset
           


cat>Job_${rootFileName}_${p}.csh<<EOF
#!/bin/tcsh
setenv X509_USER_PROXY /uscms/home/${USER}/${proxy}
source /uscmst1/prod/sw/cms/setup/cshrc prod
cd ${cmsswDir}
cmsenv
cd ${sourceDir}
/opt/d-cache/srm/bin/srmcp "file://localhost/${rootFileName}" "srm://cmssrm.fnal.gov:8443/11/store/user${destinationDir}/${rootFileName}"
echo "Copied file"
cd ${pwd}
EOF

chmod 775 $PWD/Job_${rootFileName}_${p}.csh

cat>condor_${rootFileName}_${p}<<EOF
universe = vanilla
Executable = $PWD/Job_${rootFileName}_${p}.csh
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
SSC 
Error
Log 
notify_user = schauhan@cern.ch
Queue 1
EOF

sed -e 's|SSC| Output = $(Cluster)_$(Process).stdout|'   $PWD/condor_${rootFileName}_${p} > $PWD/condor_${rootFileName}_${p}_tmp
mv $PWD/condor_${rootFileName}_${p}_tmp $PWD/condor_${rootFileName}_${p}
sed -e 's|Error| Error = $(Cluster)_$(Process).stderr|'  $PWD/condor_${rootFileName}_${p} > $PWD/condor_${rootFileName}_${p}_tmp
mv $PWD/condor_${rootFileName}_${p}_tmp $PWD/condor_${rootFileName}_${p}
sed -e 's|Log| Log = $(Cluster)_$(Process).log|'   $PWD/condor_${rootFileName}_${p} > $PWD/condor_${rootFileName}_${p}_tmp
mv $PWD/condor_${rootFileName}_${p}_tmp $PWD/condor_${rootFileName}_${p}

condor_submit condor_${rootFileName}_${p}




end

rm inputFiles.list

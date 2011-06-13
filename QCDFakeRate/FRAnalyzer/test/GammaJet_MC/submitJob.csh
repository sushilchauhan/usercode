#!/bin/tcsh
#================ Give ALL INPUTS ERE============================
setenv datafile  gjdataset.txt
setenv namefile  gjoutputfile.txt
set events = 10000                       #events per job
set tot_events = 2000000                #total events
#==================================================================

set file_length1=`wc -l $datafile | cut -c1-2`                  # Count Lines
set file_length2=`wc -l $namefile | cut -c1-2`

set p = 0

#run till the last line of the input files
while ($p != $file_length1)
@ p = ${p} + 1
set Data=`tail -n +$p ${datafile} | head -1`                  # Get name of dataset

#read the name of the output files
set DirName=`tail -n +$p ${namefile} | head -1`              # Get Name of outputfile 


#mkdir in directory in store and set permission
setenv destination_path /pnfs/cms/WAX/11/store/user/schauhan/Summer11/QCDFakeRate_PhotonTrigger/PhotonJet
setenv destination_dir ${destination_path}/${DirName}

if( ! -d ${destination_path} ) then
echo "Making directory ${destination_path}"
mkdir ${destination_path}
chmod 775 -R ${destination_path}
endif

if( ! -d ${destination_dir} ) then
echo "Making directory ${destination_dir}"
mkdir ${destination_path}/${DirName}
chmod 775 -R ${destination_path}/${DirName}
endif
#set direcotry name in store
setenv filedir /Summer11/QCDFakeRate_PhotonTrigger/PhotonJet/${DirName}

#set current directory
setenv pwd $PWD

#Give input to the crab.cfg files
cat>crab_${DirName}.cfg<<EOF
[CMSSW]
events_per_job         = ${events}
total_number_of_events = ${tot_events}
pset                   = mypatTemplate_cfg.py
datasetpath            = ${Data}
output_file            = FR_AOD_${DirName}.root

[USER]
return_data            = 0
copy_data              = 1
ui_working_dir         = ${DirName}

storage_element        = cmssrm.fnal.gov

## and the SE directory (or the mountpoint) that has to be writable from all
#### LNL SRM
storage_path           = /srm/managerv2?SFN=/11/store/user/schauhan
user_remote_dir        = ${filedir}
publish_data           = 0

[GRID]
rb                     = CERN
proxy_server           = myproxy.cern.ch
virtual_organization   = cms
retry_count            = 0

[CRAB]
scheduler              = condor
#scheduler             = glite
jobtype                = cmssw
#use_server            = 1
EOF

#Give the name of outputfile to the PATtuple files--------------------------------------------
#sed -e 's|untracked string HistOutFile=| untracked string HistOutFile="'"${DirName}.root"'"|'   ${pwd}/mypatTemplate_cfg.py > ${pwd}/mypatTemplate_cfg_tmp.py

sed -e 's|outFile          = cms.untracked.string|outFile          = cms.untracked.string('"'"FR_AOD_${DirName}.root"'"'),|'   ${pwd}/mypatTemplate_cfg.py > ${pwd}/mypatTemplate_cfg_tmp.py

cp mypatTemplate_cfg.py mypatTemplate_cfg_orig.py
mv mypatTemplate_cfg_tmp.py mypatTemplate_cfg.py
#----------------------------------------------------------------------------------------------
echo "========================================================================================================================"
echo "Submitting Job for ${Data} with Ouput file Name ${DirName}"
echo "========================================================================================================================"

#sleep for 10 minuts

#create and submit the job---------------------------------------------
if(-f log_Submit_${DirName})then
crab -create -submit -cfg crab_${DirName}.cfg >> log_Submit_${DirName}
else
crab -create -submit -cfg crab_${DirName}.cfg > log_Submit_${DirName}
endif
#--------------------------------------------------------------------

#now put back the original file at its place
cp mypatTemplate_cfg.py junk/mypatTemplate_${DirName}_cfg.py
mv mypatTemplate_cfg_orig.py mypatTemplate_cfg.py 
cp crab_${DirName}.cfg junk/crab_${DirName}.cfg
#---------------------------------------------------------------------
end




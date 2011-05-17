#!/bin/tcsh
#--------------------Change for your user account---------
#setenv X509_USER_PROXY /uscms/home/sandhya/x509up_u45353
#---------------------------------------------------------
setenv pwd $PWD
echo ${pwd}


#----------SET THESE VARIABLE HERE----------------------------------------
setenv sourceDir /uscms_data/d2/sushil/CMSSW/MonoPhoton/CMSSW_3_8_7/src
set NFiles = 100
setenv LHE_FILE z1nunugamma_sm_nlo       # Name without .lhe extension
setenv output_location /uscmst1b_scratch/lpc1/3DayLifetime/sushil
#------------------------------------------------------------------------- 



set p = 0

#Run Below for N files
while ($p != $NFiles)
@ p = ${p} + 1


cat>Pythia8HadAndSim_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_${p}.py<<EOF
import FWCore.ParameterSet.Config as cms                                                                                                                                           
from Configuration.Generator.PythiaUESettings_cfi import *
 
process = cms.Process('HLT')
 
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.VtxSmearedRealistic7TeVCollision_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')
 
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.222 $'),
    #annotation = cms.untracked.string('Pythia8HadAndSim_cfi.py nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(
 
)
# Input source
process.source = cms.Source("LHESource",
    fileNames = cms.untracked.vstring('file:/uscms_data/d2/carley/photon/CMSSW_3_8_4/src/BaurLHE/LHE/LO/LHC_FILE_${p}.lhe')
)
 
# Output definition
 
process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('file:${output_location}/znulo1_sm_step1_${p}.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)
 
# Additional output definition
 
# Other statements
process.GlobalTag.globaltag = 'START38_V13::All'
process.generator  = cms.EDFilter("Pythia8HadronizerFilter",
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    comEnergy = cms.double(7000.),
    #LHEInputFileName = cms.untracked.string('ttbar.lhe'),
    PythiaParameters = cms.PSet(
        pythia8_dhidas = cms.vstring(''),
        parameterSets = cms.vstring('pythia8_dhidas')
    ),
 
    # dhidas
    jetMatching = cms.untracked.PSet(
      scheme = cms.string("Pythia8BaurWGammaNLO"),
      applyMatching = cms.bool(True)
#I think the gen pT jet cut should be defined in here
    )
 
)
 
 
# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.Path(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)
 
# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.RAWSIMoutput_step])
# special treatment in case of production filter sequence  
for path in process.paths:
    getattr(process,path)._seq = process.generator*getattr(process,path)._seq         

EOF






cat>Job_${p}.csh<<EOF
#!/bin/tcsh
source /uscmst1/prod/sw/cms/setup/cshrc prod
cd $sourceDir
cmsenv
cd ${pwd}
cmsRun ${pwd}/Pythia8HadAndSim_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_${p}.py
EOF





chmod 775 ${pwd}/Job_${p}.csh





cat>condor_${p}<<EOF
universe = vanilla
Executable = ${pwd}/Job_${p}.csh
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
SSC 
Error
Log 
notify_user = sushil@fnal.gov
Queue 1
EOF

sed -e 's|SSC| Output = $(Cluster)_$(Process).stdout|'   ${pwd}/condor_${p} > ${pwd}/condor_${p}_tmp
mv ${pwd}/condor_${p}_tmp ${pwd}/condor_${p}
sed -e 's|Error| Error = $(Cluster)_$(Process).stderr|'  ${pwd}/condor_${p} > ${pwd}/condor_${p}_tmp
mv ${pwd}/condor_${p}_tmp ${pwd}/condor_${p}
sed -e 's|Log| Log = $(Cluster)_$(Process).log|'   ${pwd}/condor_${p} > ${pwd}/condor_${p}_tmp
mv ${pwd}/condor_${p}_tmp ${pwd}/condor_${p}

#condor_submit condor_${p}



end


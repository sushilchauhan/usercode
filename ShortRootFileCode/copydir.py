import commands, string, re

get = commands.getoutput('lcg-ls -b -D srmv2 --verbose "srm://f-dpm001.grid.sinica.edu.tw:8446/srm/managerv2?SFN=/dpm/grid.sinica.edu.tw/home/cms/store/user/cmkuo/pythia_znunug_413_GEN_SIMv1/pythia_znunug_413_RECO_v1/53487bfe0f4ffcbcde83ffeb864c5dec/" | grep 53487bfe0f4ffcbcde83ffeb864c5dec/.*$')

lines = string.split(get, "\n")
for line in lines:
	print line
	# Search the file name for the string 'askew' -- this is a basic way of selecting certain files
	if re.search('STEP2', line):
		parts = string.split(line, "/")
		#Change below to point to the appropriate directory on the storage element
		command = 'lcg-cp --verbose -b -D srmv2 "srm://f-dpm001.grid.sinica.edu.tw:8446/srm/managerv2?SFN=' + line + '" "srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/11/store/user/sushil/Summer11/MC/RECO/ZNuNu_Gamma/' + parts[-1] + '"'
		print command
		#Execute the above command (comment out the line below to test)
		output = commands.getoutput(command)

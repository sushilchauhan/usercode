#!/bin/tcsh                                                                                                                                                                            
set file_length=${2}
 
set p = 0
set t = ${1}
 
#run till the last line of the input files
while ($p != $file_length)
@ p = ${p} + 1
condor_rm ${t}.${p}
echo "killed " ${t}.${p}
end


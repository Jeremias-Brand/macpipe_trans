#!/bin/bash

echo "These are the versions of software used for run: {params.time}" 
date                                                                  
echo "Operating system:"                                              
cat /etc/centos-release                                              
echo "linux kernel:"                                                
uname -a                                                              
echo "gcc version (careful some software might be using different libs)"og}
gcc --version | grep gcc                                              
snap-aligner 2>&1 | grep version                                      
salmon -version                                                       
R --version | grep "R version"                                        
python --version                                                      
echo "software installed via anaconda:"                               
conda list                                                            
echo "python modules:"                                                
pip freeze                                                            
hmmscan -h | grep -E '^\#\ H'                                         
blastn -h | grep Nucleotide                                           
embossversion                                                         
echo "perl version: "                                                 
perl -e 'print $];'                                                   
echo "BUSCO version: "                                                
$3 --version                                              
echo "transrate version:" 
$2 --version                
echo $1                                                 
echo "HARDCODED"                                                      
echo "TransDecoder   2.0.1"                                           
echo "Trimmomatic 0.36"                                               


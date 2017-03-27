# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 15:06:17 2016

@author: Jeremias Brand

Takes assemblies files from transrate run as input and concatenates them.
"""

assemblies = ["Mcla2_subsample.assemblies.csv", "Mcla_subsample.assemblies.csv"]
#generate first line we need in outfile:
with open("cat.assemblies", "w") as outfile:
    with open(assemblies[0], "r") as infile:
        lines = infile.readlines()
        outfile.write(lines[0])
    
# now we open the outfile in appending mode to add all the results
with open("cat.assemblies", "a") as outfile:
    for assembly in assemblies:
        with open(assembly, "r") as infile:
            lines = infile.readlines()
            outfile.write(lines[1])




#with open("Mcla2_subsample.assemblies.csv", "r") as infile:
#    lines = infile.readlines()
#    print(lines)
#    print("second line:")
#    print(lines[1])
#    
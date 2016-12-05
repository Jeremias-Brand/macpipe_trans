#! /home/jeremias/.linuxbrew/bin/R

print("this is object one:")
print(snakemake@input[[1]])
file.create("r.sentinel")

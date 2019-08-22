library(scImpute)
args <- commandArgs(TRUE)
scimpute(count_path=args[1], infile='csv', outfile='csv',out_dir = "./",drop_thre = 0.5, ncores = 10, Kcluster=strtoi(args[2]))
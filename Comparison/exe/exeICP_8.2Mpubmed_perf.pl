#!/bin/perl

$ALGO = "ICP"; $algo = "icp";
$HOMEDIR = "../../";
$CompDIR = "$HOMEDIR/Comparison";
$data_type = "pubmed";
$data_size = "8_2M";
$data_name = "$data_type.$data_size.db";
$DATADIR = "$HOMEDIR/dataset";
	if(! -d $DATADIR){printf("No DATADIR exists: %s\n",$DATADIR);}
$BINDIR = "$CompDIR/bin";
	if(! -d $BINDIR){printf("No BINDIR exists: %s\n",$BINDIR);}
$OUTDIR = "$CompDIR/Log";

$perf = "perf stat -e cpu-cycles -e instructions -e branches -e branch-misses -e L1-dcache-loads -e L1-dcache-load-misses -e LLC-loads -e LLC-load-misses -e LLC-stores -e LLC-store-misses -e cache-references -e cache-misses ";

$NumThreads = 50;

$bin = "$BINDIR/$algo";
	if(! -x $bin){printf("No bin file: %s\n",$bin);}

$input = "$DATADIR/$data_name";
	if(! -r $input){printf("No input exists.\n");exit(1);}

$seed = 0;
$seedN = $seed; $seedK = 100+$seed;

#@K = (80000,40000,20000,10000);
@K = (80000);
foreach $num_means (@K){
$fname = "$algo.N$data_size.K$num_means.s$seed";
$out_assign = "$OUTDIR/$fname.assign";
$out_mem = "$OUTDIR/$fname.mem_info";
$out_log = "$OUTDIR/$fname.log";

system("$perf $bin $NumThreads $seedN $seedK $num_means $input $out_assign $out_mem 2> $out_log");

}

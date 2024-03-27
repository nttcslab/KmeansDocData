#!/bin/perl

$ALGO = "AF-ICP"; $algo = "itr5_aficp";
$HOMEDIR = "../";
$data_type = "pubmed";
$data_size = "8.2M";
$data_name = join("_",$data_size,$data_type);
$DATADIR = "$HOMEDIR/dataset";
	if(! -d $DATADIR){printf("No DATADIR exists: %s\n",$DATADIR);}
$BINDIR = "$HOMEDIR/bin";
	if(! -d $BINDIR){printf("No BINDIR exists: %s\n",$BINDIR);}
$OUTDIR = "$HOMEDIR/Log";

$perf = "perf stat -e cpu-cycles -e instructions -e branches -e branch-misses -e L1-dcache-loads -e L1-dcache-load-misses -e LLC-loads -e LLC-load-misses -e LLC-stores -e LLC-store-misses -e cache-references -e cache-misses ";

$NumThreads = 50;

$bin = "$BINDIR/$algo";
	if(! -x $bin){printf("No bin file: %s\n",$bin);}

$ThTermMin = 122000;
$thterm = $ThTermMin; $thterm =~ s/000$/K/;
$ThValMin = 0.020; $ThValMax = 0.044;
$ThValStep1 = 0.002; $ThValStep2 = 0.001;

$input = "$DATADIR/$data_name.db";
	if(! -r $input){printf("No input exists.\n");exit(1);}

$seed = 0;
$seedN = $seed; $seedK = 100+$seed;

@K = (10000);
foreach $num_means (@K){
$fname = "$algo.N$data_size.K$num_means.s$seed";
$out_assign = "$OUTDIR/$fname.assign";
$out_mem = "$OUTDIR/$fname.mem_info";
$out_log = "$OUTDIR/$fname.log";

system("$perf $bin $NumThreads $ThTermMin $ThValMin $ThValMax $ThValStep1 $ThValStep2 $seedN $seedK $num_means $input $out_assign $out_mem 2> $out_log");

}

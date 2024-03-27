#!/bin/perl

$DATADIR = "./";
	if(! -d $DATADIR){printf("No DATADIR: %s\n",$DATADIR);exit(1);}
$fname = "docword.pubmed.txt";
$inf = "$DATADIR/$fname";
	if(! -r $inf){printf("No input:%s\n",$inf);exit(1);}
$outf = "$DATADIR/8.2M_pubmed.db";

undef(%TF); undef(%NUM);
open(IN,$inf);
	$N = <IN>; chop($N); $D = <IN>; chop($D);
	$nnz = <IN>; chop($nnz);
	$num_words = 0; $prev_docID = 1; $max_wordID = 0;
	while(<IN>){
		chomp; undef(@tmp);
		@tmp = split(/\s+/,$_);
		$docID = $tmp[0];
		$wordID = $tmp[1];
		if($wordID > $max_wordID){$max_wordID = $wordID;}
		$tf = $tmp[$#tmp];
		$TF{$docID}{$wordID} = $tf;
		if($docID == $prev_docID){
			$num_words++;
			if(eof(IN)){
				$NUM{$docID} = $num_words;
			}
		}else{
			$NUM{$docID-1} = $num_words;
			$num_words = 1;	
			$prev_docID = $docID
		}
	}
close(IN);
$N = $docID;
$D = $max_wordID;

open(OUT,">$outf");
printf(OUT "%d %d 1\n",$N,$D);
foreach $docID (sort{$a<=>$b} keys %TF){
	printf(OUT "%d",$NUM{$docID});
	foreach $wordID (sort{$a<=>$b} keys %{$TF{$docID}}){
		printf(OUT " %d:%d",$wordID,$TF{$docID}{$wordID});
	}
	printf(OUT " 1 1\n");
}
close(OUT);

#!/usr/bin/perl/ -w
$|++; #---turn on the auto flush for the progress bar
use strict;
use Math::Random; #---For generating normal distribution random numbers
use File::Path; #---for removing tmp GNUPLOT dat files
use List::Util 'shuffle';#--- for shuffling an array
use POSIX qw(log10);
use Time::HiRes qw( time );

######################################################################################################################################################
#
#	Description
#		This is a perl script to generate random reads from a genomic fasta file. Both single and pair end reads are supported. 
#
#	Input
#		
#		--fasta=		filename; name of the fasta file;
#		--readNum=		integer; read number in millions, default 10;
#		--readLen=		integer; read length, default = 50;
#		--SPEnd=		single or pair; single end or pair end library;
#		--fragLenMean=	integer; the mean of the normal distribution of fragment length, note that mean +- 2 SD includes ~95% of the values;
#		--fragLenSD=	integer; the standard devation of the normal distribution of fragment length, note that mean +- 2 SD includes ~95% of the values;
#		--maxFragLen=	integer; enforce an upper limit for the normal distribution of fragment length;
#		--minFragLen=	integer; enforce a lower limit for the normal distribution of fragment length;
#		--sam=	yes or no; to output a sam alignment file that specify the location of the simulated read; The purpose of doing this is to visualize whether some reads are missed during the alignment process, e.g. reads derived highly repetitive regions;
#
#	Output
#		1. a fragment length distribution plot. Tha fragment length with be drawn from this distribution;
#		2. a fasta file (or two fasta files in the case of pair-end. p.s. the naming of the pairend header might have problems, will look into it later on);
#
#	Usage
#		perl gReadGenerator_v0.2.pl --fasta=clean.EhistolyticaGenomic_AmoebaDB-1.0WithrDNA.fasta --readNum=10 --readLen=20 --SPEnd=single --fragLenMean=175 --fragLenSD=15 --maxFragLen=200 --minFragLen=150 --sam=yes
#
#
#	Assumptions
#
#	Version history
#
#		v0.1:	-debut;
#		
#		v0.2: -the option --sam= is added. It outputs a sam file which contains the alignment of the  
#
######################################################################################################################################################

#==========================================================Main body starts==========================================================================#

#1----------Read the parameters----------#
use vars qw ($fasta $readNum $readLen $SPEnd $fragLenMean $fragLenSD $maxFragLen $minFragLen $sam $paraTag);
($fasta, $readNum, $readLen, $SPEnd, $fragLenMean, $fragLenSD, $maxFragLen, $minFragLen, $sam, $paraTag) = readParameters();

#2----------check the parameters----------#
checkParameters();

printCMDLogOrFinishMessage("CMDLog");

#3----------read the fasta----------#
my $fastaHsh_ref = readMultiFastaFile();

#4----------generate fragment length distribution----------#
my $fragLenDstbtnAry_ref = generateFragLenDistributionArray();

#5----------generate the read---------#
generateReads($fastaHsh_ref, $fragLenDstbtnAry_ref);

printCMDLogOrFinishMessage("finishMessage");
exit;
#========================================================= Main body ends ===========================================================================#

########################################################################## readParameters
sub readParameters {
	
	$readNum = 10;
	$readLen = 50;
	$SPEnd = "single";
	$fragLenMean = 175;
	$fragLenSD = 50;
	$maxFragLen = 200;
	$minFragLen = 150;
	$sam = "no";
	
	foreach my $param (@ARGV) {
		if ($param =~ m/--fasta=/) {$fasta = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--readNum=/) {$readNum = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--readLen=/) {$readLen = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--SPEnd=/) {$SPEnd = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--fragLenMean=/) {$fragLenMean = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--fragLenSD=/) {$fragLenSD = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--maxFragLen=/) {$maxFragLen = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--minFragLen=/) {$minFragLen = substr ($param, index ($param, "=")+1);}
		elsif ($sam =~ m/--sam=/) {$sam = substr ($param, index ($param, "=")+1);}
	}
	
	my $fastaxExt = $fasta;
	$fastaxExt =~ s/\.\w+$//; 
	my $paraTag = $fastaxExt."_rNum.".$readNum."_rLen.".$readLen."_".$SPEnd."_frag.".$fragLenMean.".".$fragLenSD.".".$maxFragLen.".".$minFragLen;

	return ($fasta, $readNum, $readLen, $SPEnd, $fragLenMean, $fragLenSD, $maxFragLen, $minFragLen, $sam, $paraTag);
}
########################################################################## checkParameters
sub checkParameters {

	die "The read number is too small. Better be something > 0.1 million.\n" if ($readNum < 0.1);
	die "The read number is too large. Better be something < 1000 million.\n" if ($readNum > 1000);
	die "The read length is too short. Better be something > 20 nt.\n" if ($readLen < 20);
	die "The read length is too long. Better be something < 500 nt.\n" if ($readLen > 500);
	die "The SPEnd options has to be \"single\" or \"pair\".\n" if (($SPEnd ne "single") and ($SPEnd ne "pair"));
	die "minFragLen > fragLenMean or maxFragLen < fragLenMean, which is impossible. Please reset.\n" if (($minFragLen > $fragLenMean) or ($maxFragLen < $fragLenMean));
	die "fragLenSD has to be smaller than 1/3 of fragLenMean, please reset.\n" if ($fragLenSD > ($fragLenMean*0.3));
	die "readLen is longer than fragLenMean, please reset.\n" if ($readLen > $fragLenMean);
	
}
########################################################################## readMultiFastaFile
sub readMultiFastaFile {

	my ($seq, $seqName, %fastaHsh);
	my $i = 0;
	print "Reading $fasta into a hash.\n";	
	open (INFILE, $fasta);
	chomp (my $theCurntLine = <INFILE>); #get the first line
	while (my $theNextLine = <INFILE>) {
		chomp $theNextLine;
		
		#---Only two types of line in current line, the header or seq
		if ($theCurntLine =~ m/^>/) {#-- header line
			my @theLineSplt = split (/\|/, $theCurntLine);
			$seqName = $theLineSplt[0]; #---get the second tag
			$seqName =~ s/ //g; #---remove space
			$seqName =~ s/>//g; #---remove space
		} else {#--seq line
			$seq = $seq.$theCurntLine;
		}
		
		#---check if next line has a > or that's the end of file
		if ($theNextLine =~ m/^>/) {
			$fastaHsh{$seqName} = $seq;
			$seq = "";
		} elsif (eof(INFILE)) {#---this is the last line
			$seq = $seq.$theNextLine;
			$fastaHsh{$seqName} = $seq;
		}
		
		#---next line becomes current line
		$theCurntLine = $theNextLine;
	}
	close INFILE;

	return (\%fastaHsh);
}
########################################################################## generateFragLenDistributionArray
sub generateFragLenDistributionArray {

	print "\nGenerating fragment length distribution by sampling 1 million random numbers.\n";

	#----random 1 million numbers from the desired range
	my $fragLenDstbtnAry_ref = randNonZeroNumNormDistWithCuttOffs($fragLenMean, $fragLenSD, 1000000, $maxFragLen, $minFragLen);
	
	my @fragLenDstbtnAry = @{$fragLenDstbtnAry_ref};
	
	#---store the percentage of frag length into a hash and plot a XY plot
	my %fragLenProbHsh;
	my $totalFragLenNum = @fragLenDstbtnAry;
	foreach my $fragLen (@fragLenDstbtnAry) {
		$fragLenProbHsh{$fragLen}++;
	}

	foreach my $fragLen (keys %fragLenProbHsh) {
		$fragLenProbHsh{$fragLen} = $fragLenProbHsh{$fragLen}/$totalFragLenNum;
	}
	
	my $fragLenProbHsh_ref = \%fragLenProbHsh;
	
	GNUPLOTHshXYScatter($fragLenProbHsh_ref, "$paraTag.fragLenDstbtn.pdf", "fragment length (nt)", "probability", "linear", "linear", "fragment length distribution from 10000 random numbers");
	
	return (\@fragLenDstbtnAry);
	
}
########################################################################## randNonZeroNumNormDistWithCuttOffs
sub randNonZeroNumNormDistWithCuttOffs {

	my $mean = $_[0];
	my $sd = $_[1];
	my $num = $_[2]; #---number of random numbers needed
	my $upperLimit = $_[3];
	my $lowerLimit = $_[4];

	my $procNum = 0;
	my $interval = int ($num/50);
	my $printInterval = 0;
	
	print "|--------------------------------------------------|100%\n";
	print "|";
	
	my @normRandNumAry;
	
	for (1..$num) {
		my $tmpRandNum = -1;
		while (($tmpRandNum <= 0) or ($tmpRandNum >= $upperLimit) or ($tmpRandNum <= $lowerLimit)) {
			$tmpRandNum = int (random_normal(1, $mean, $sd));
		}
		$procNum++; 
		if (($procNum == $interval) and ($printInterval < 50)) {
			print "-";
			$procNum = 0;
			$printInterval++;
		}
		push (@normRandNumAry, $tmpRandNum);
	}
	print "|done\n\n";
	
	return (\@normRandNumAry);
}
########################################################################## GNUPLOTHshXYScatter
sub GNUPLOTHshXYScatter {

	my %XYHsh = %{$_[0]};
	my $filePath = $_[1];
	my $xlable = $_[2];
	my $ylable = $_[3];
	my $xscale = $_[4];
	my $yscale = $_[5];
	my $title = $_[6];
	
	
	my $GNULogXCmd = "";
	$GNULogXCmd = "set logscale x" if ($xscale eq "log");
	my $GNULogYCmd = "";
	$GNULogYCmd = "set logscale y" if ($yscale eq "log");

	print "Running GNUPLOTHshXYScatter for $filePath.\n";
	
	#---creat a tmp file
	open (TMPFILE, ">tmp.dat");
	for my $x (sort {$a <=> $b} keys %XYHsh) {
		print TMPFILE $x."\t".$XYHsh{$x}."\n";
	}
	close TMPFILE;
	
	#---do the GNUPLOT
	open (GNUPLOT, "|gnuplot");
	print GNUPLOT <<EOPLOT;
	set terminal postscript color solid
	set output "| ps2pdf - $filePath 2>/dev/null";
	unset logscale x; 
	unset logscale y; 
	$GNULogXCmd;
	$GNULogYCmd;
	set xlabel "$xlable";
	set ylabel "$ylable";
	set title "$title";
	set nokey;
   	plot 'tmp.dat' using 1:2 with lines;
EOPLOT
	close(GNUPLOT);
	rmtree(['tmp.dat'], 0, 1); #---non-verbose removal of tmp file
}
########################################################################## generateReads
sub generateReads {
	
	my %fastaHsh = %{$_[0]};
	my @fragLenDstbtnAry = @{$_[1]};
	
	#----assign the proportion of reads to each of the cntg
	my $totalSeqLen = 0;
	my %readNumForEachCntg;
	
	foreach my $seqName (keys %fastaHsh) {#---get the total length
		$totalSeqLen = $totalSeqLen + (length $fastaHsh{$seqName});
	}
	
	foreach my $seqName (sort {$a cmp $b} keys %fastaHsh) {
		$readNumForEachCntg{$seqName} = int (((length $fastaHsh{$seqName}) / $totalSeqLen)*($readNum*1000000));
		print "$readNumForEachCntg{$seqName} reads will be generated for contig $seqName.\n";
		die "contig $seqName is having zero reads assigned. Please increase the read number.\n" if ($readNumForEachCntg{$seqName} <= 0);
	}
	
	open (FASTAOUT1, ">$paraTag.1.fasta");
	open (FASTAOUT2, ">$paraTag.2.fasta") if ($SPEnd eq "pair");
	
	#---open the sam file and print the sam header
	my $qualityString = "";
	for my $i ((1..$readLen)) {#----quality string for SAM
		$qualityString .= "I";
	}

	if ($sam eq "yes") {
		open (SAMOUT, ">$paraTag.original.sam");
		print SAMOUT "\@HD\tVN:1.0\tSO:unsorted\n";
		foreach my $seqName (sort {$a cmp $b} keys %fastaHsh) {
			my $cntgLen = length $fastaHsh{$seqName};
			print SAMOUT "\@SQ\tSN:$seqName\tLN:$cntgLen\n"
		}
		my $command = join " ", @ARGV;
		print SAMOUT "\@PG\tID:gReadGenerator\tVN:0.2\tCL:perl gReadGenerator_v0.2.pl $command\n";
	}

	my $readProc = 0;
	my $progCount = 0;
	my $totalTimeSpent = 0;
	my $averageTimeSpent;
	print "Start generating the fragments.\n";
	my ($intervalStart, $intervalEnd);
	$intervalStart = time();

	#---Looping through each contig
	foreach my $seqName (sort {$a cmp $b} keys %fastaHsh) {

		my $cntgSeq = $fastaHsh{$seqName};
		my $cntgLen = length $cntgSeq;

		#---loop through each fragment
		for my $i (1..$readNumForEachCntg{$seqName}) {
		
			$readProc++;
			$progCount++;
			
			if ($progCount == 100000) {

				print "$readProc fragments have been processed. "; 
				$progCount=0;
				$intervalEnd = time();
				my $timeElapsed = $intervalEnd - $intervalStart;
				$timeElapsed = sprintf ("%.2f", $timeElapsed);
				$totalTimeSpent += $timeElapsed;
				$averageTimeSpent = sprintf ("%.2f", ($totalTimeSpent/$readProc)*100000);
				
				my $estimatedEnd = ((($readNum*1000000)- $readProc)*$averageTimeSpent)/100000;
				$estimatedEnd = sprintf ("%.2f", $estimatedEnd/60);
				print "Average Time to process 100000 fragment is ".$averageTimeSpent." sec. Estimated to end in ".$estimatedEnd." mins.\n";
				$intervalStart = time();
			}
			
			my $fragLen = $fragLenDstbtnAry[int (rand @fragLenDstbtnAry)];
			my $fragLenHalf = int ($fragLen/2);
			my ($fragMidPos, $fragStartPos, $fragEndPos);

			my $validPos = "no";
			my $fragSeqPlus;

			while ($validPos eq "no") {#---loop until not out of ctng rng
				my $outCngtRng = "yes";
				my $seqContainNs = "yes";
				$fragMidPos = int (rand $cntgLen); #---rand from 0 to $cntgLen-1;
				$fragStartPos = $fragMidPos - $fragLenHalf;
				$fragEndPos = $fragStartPos + $fragLen;
				$outCngtRng = "no" if (($fragEndPos < $cntgLen) and ($fragStartPos >= 0));
				$fragSeqPlus = substr $cntgSeq, $fragStartPos, $fragLen;
			 	$seqContainNs = "no" if ($fragSeqPlus !~ m/[^ATGCatgc]/);
			 	$validPos = "yes" if (($outCngtRng eq "no") and ($seqContainNs eq "no"));
			}
			
			#---get the fragseq and reverse complement the seq
			my $fragSeqMinus = reverse $fragSeqPlus;
			$fragSeqMinus =~ tr/ACGTacgt/TGCAtgca/;
			
			#---rand for the end in which the adaptor will be ligated
			my $dir = "+";
			$dir = "-" if (rand() >= 0.5);
			
			#----get the read seq from the end of the fragment
			my ($readSeq1, $readSeq2, $readStartPos1, $readStartPos2, $dir1, $dir2, $read1Flag, $read2Flag, $read1FragLen, $read2FragLen, $samSeq1, $samSeq2);
			
			if ($dir eq "+") {
				
				$dir1 = "+";
				$dir2 = "-";
				$readSeq1 = substr $fragSeqPlus, 0, $readLen;
				$readSeq2 = substr $fragSeqMinus, 0, $readLen;

				$samSeq1 = $readSeq1;
				$samSeq2 = reverse $readSeq2;
				$samSeq2 =~ tr/ACGTacgt/TGCAtgca/;

				$readStartPos1 = $fragStartPos + 1;
				$readStartPos2 = $fragEndPos - $readLen + 1;
				$read1Flag = 0;
				
				if ($SPEnd eq "pair") {
					$read1Flag = 1+  2+  0+  0+  0+   32+  64+  0; #= 99;
					$read2Flag = 1+  2+  0+  0+  16+   0+   0+  128; #= 147
					$read1FragLen = $fragEndPos - $fragStartPos;
					$read2FragLen = $fragStartPos - $fragEndPos;
				}
				
			} elsif ($dir eq "-") {

				$dir1 = "-";
				$dir2 = "+";
				$readSeq1 = substr $fragSeqMinus, 0, $readLen;
				$readSeq2 = substr $fragSeqPlus, 0, $readLen;

				$samSeq1 = reverse $readSeq1;
				$samSeq1 =~ tr/ACGTacgt/TGCAtgca/;
				$samSeq2 = $readSeq2;

				$readStartPos1 = $fragEndPos - $readLen + 1;
				$readStartPos2 = $fragStartPos + 1;
 				$read1Flag = 16;

				if ($SPEnd eq "pair") {
					$read1Flag = 1+  2+  0+  0+   0+  32+   0+  128; #= 163
					$read2Flag = 1+  2+  0+  0+  16+   0+  64+  0; #= 83
					$read1FragLen = $fragStartPos - $fragEndPos;
					$read2FragLen = $fragEndPos - $fragStartPos;
				}
			}
		
#			0x0001 1 paired end is true
#			0x0002 2 both mates mapped  
#			0x0004 4 the query sequence itself is unmapped 
#			0x0008 8 the mate is unmapped  
#			0x0010 16 strand of the query (0 for forward; 1 for reverse strand) 
#			0x0020 32 strand of the mate  (0 for forward; 1 for reverse strand)
#			0x0040 64 the read is the ﬁrst read in a pair  
#			0x0080 128 the read is the second read in a pair 
#			0x0100 256 the alignment is not primary (a read having split hits may have multiple primary alignment records) 
			
			#---output the sequence
			my $readName1 = $seqName."_".$readLen."M"."_".$dir1.$readStartPos1."_".$dir2.$readStartPos2."/1";
			my $readName2 = $seqName."_".$readLen."M"."_".$dir1.$readStartPos1."_".$dir2.$readStartPos2."/2";
			
			print FASTAOUT1 "@".$readName1."\n";
			print FASTAOUT1 $readSeq1."\n";
			print FASTAOUT1 "+\n";
			print FASTAOUT1 $qualityString."\n";
			
			my @samLineAry = ();
			my $samLineString;

			if (($sam eq "yes") and ($SPEnd ne "pair")) {
				#Print sam file
				#DS571145_50M_+77324_-77454	99	DS571145	77324	255	50M	=	77454	180	TGTGGGATAATGATGAAGATGAACCTATTAAAGATGATTTGAGTGATGGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	XA:i:0	MD:Z:50	NM:i:0
		
				@samLineAry = ();
				$samLineAry[0] = $readName1;
				$samLineAry[1] = $read1Flag;
				$samLineAry[2] = $seqName;
				$samLineAry[3] = $readStartPos1;
				$samLineAry[4] = 255;
				$samLineAry[5] = $readLen."M";
				$samLineAry[6] = "*";
				$samLineAry[7] = "*";
				$samLineAry[8] = "*";
				$samLineAry[9] = $samSeq1;
				$samLineAry[10] = $qualityString;
				$samLineAry[11] = "XA:i:0";
				$samLineAry[12] = "MD:Z:50";
				$samLineAry[13] = "NM:i:0";
				$samLineString = join "\t", @samLineAry;
				print SAMOUT $samLineString."\n";
			}

			if ($SPEnd eq "pair") {
				print FASTAOUT2 "@".$readName2."\n";
				print FASTAOUT2 $readSeq2."\n";
				print FASTAOUT2 "+\n";
				print FASTAOUT2 $qualityString."\n";
				if ($sam eq "yes") {
					#Print sam file
					#DS571145_50M_+77324_-77454	99	DS571145	77324	255	50M	=	77454	180	TGTGGGATAATGATGAAGATGAACCTATTAAAGATGATTTGAGTGATGGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	XA:i:0	MD:Z:50	NM:i:0

					my $readName1Out = $readName1;
					$readName1Out =~ s/\/1//;
					@samLineAry = ();
					$samLineAry[0] = $readName1Out;
					$samLineAry[1] = $read1Flag;
					$samLineAry[2] = $seqName;
					$samLineAry[3] = $readStartPos1;
					$samLineAry[4] = 255;
					$samLineAry[5] = $readLen."M";
					$samLineAry[6] = "=";
					$samLineAry[7] = $readStartPos2;
					$samLineAry[8] = $read1FragLen;
					$samLineAry[9] = $samSeq1;
					$samLineAry[10] = $qualityString;
					$samLineAry[11] = "XA:i:0";
					$samLineAry[12] = "MD:Z:50";
					$samLineAry[13] = "NM:i:0";
					$samLineString = join "\t", @samLineAry;
					print SAMOUT $samLineString."\n";

					my $readName2Out = $readName2;
					$readName2Out =~ s/\/2//;
					@samLineAry = ();
					$samLineAry[0] = $readName2Out;
					$samLineAry[1] = $read2Flag;
					$samLineAry[2] = $seqName;
					$samLineAry[3] = $readStartPos2;
					$samLineAry[4] = 255;
					$samLineAry[5] = $readLen."M";
					$samLineAry[6] = "=";
					$samLineAry[7] = $readStartPos1;
					$samLineAry[8] = $read2FragLen;
					$samLineAry[9] = $samSeq2;
					$samLineAry[10] = $qualityString;
					$samLineAry[11] = "XA:i:0";
					$samLineAry[12] = "MD:Z:50";
					$samLineAry[13] = "NM:i:0";
					$samLineString = join "\t", @samLineAry;
					print SAMOUT $samLineString."\n";

				}#---end of if ($sam eq "yes") {
			}#---end of if ($SPEnd eq "pair") {
		}#---end of for my $i (1..$readNumForEachCntg{$seqName}) {
	}#---end of foreach my $seqName (sort {$a cmp $b} keys %fastaHsh) {	
}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptNameXext = $0;
		$scriptNameXext =~ s/\.\w+$//;
		open (CMDLOG, ">>$scriptNameXext.cmd.log.txt"); #---append the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}
	
}

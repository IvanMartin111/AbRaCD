#!/usr/bin/perl
# Refine RCD-generated loops and summarize their Rosetta scores.

use strict;
use File::Basename;



# Helper scripts are resolved from the directory where this script is run.
my $minimizer = "loops_energy13.py";
my $relax = "monrelax4.py";
my $insert = "insertpdb.pl";
my $inpref = "../";
my $col_rosettaE = 9; # Column index (for -k option of sort command)
my $col_rmsd = 12; # Column index of PyRosetta refined RMSD
my $run = 1; # Controls helper commands used for insertion and output-file processing.

# PyRosetta settings for native loops and modeling workflows.
my $optmin = "--nmin 5 --initneighbor -1 --minneighbor 0 --neighbor -1 --irepack 0 --mintype dfpmin_armijo_nonmonotone --nrep 1 --fas=0.05,0.30,0.60,1.00 --ntpsi --verbose error ";
my $optrelax = "--nmin 5 --initneighbor 10 --nminneigh 100 --irepack 0 --mintype dfpmin_armijo_nonmonotone --fas=0.05,0.30,0.60,1.00 --ntpsi --deleteloop ";
my $optrelax_withloop = "--nmin 5 --initneighbor 10 --nminneigh 100 --irepack 0 --mintype dfpmin_armijo_nonmonotone --ntpsi";
my $optmin2 = "--nmin 5 --initneighbor -1 --minneighbor 5 --neighbor -1 --irepack 0 --mintype dfpmin_armijo_nonmonotone --nrep 1 --fas=0.05,0.30,0.60,1.00 --ntpsi --verbose error ";


# Validate the required input arguments.
if( !($#ARGV == 5 || $#ARGV == 4) )
{
	print "USAGE:\n\t$0 <RCD_input_file> [native|modeling] [outsuff] [np] [\"parameters\"]>(Optional)\n\n";
	print "\t<RCD_input_file>  --> RCD's input text file indicating: loop PDB, start and end indices, chain, and sequence.\n";
	print "\t[native|modeling] --> Choose between \"native\" or \"modeling\" scenarios.\n";
	print "\t[outsuff]         --> Suffix for output files. E.g: \"_ref2015\"\n";
	print "\t[np]              --> Number of processes for MPI parallelization.\n";
	print "\t<parameters>      --> (Optional) Additional parameters for PyRosetta script. E.g: \"--energy ref2015\"\n";
	print "\n\tNOTE that current PyRosetta parameters are: \n";
	print "\t*optmin ==> $optmin\n";
	print "\t*optrelax ==> $optrelax\n";
	print "\t*optmin2 ==> $optmin2\n";
	print "\nDESCRIPTION:\n\tRCD+ loop refinement for \"native\" , \"modeling\" ,\"crista2model\" or \"stats\" scenarios.\n";
	print "\nWARNING:\n\tEnsure that the scripts: monrelax4.py loops_energy13.py monfuncs.py insertpdb.pl are on PATH.\n";
	print "\nEXAMPLE:\n";
	print "\ttime refine_loop.pl ../mini8.txt native _ref2015nat 4 \"--energy ref2015\"\n";
	print "\ttime refine_loop.pl ../mini8.txt modeling _ref2015mod 4 \"--energy ref2015\"\n";
	print "\ttime refine_loop.pl ../mini8.txt stats _ref2015mod 4 \"--energy ref2015\"\t#scenario for only recalculate thet stats files without using rosetta\n";
	print "\nCRISTALTOMODEL steps:\n";
	print "\t1-	time refine_loop.pl ../mini8.txt crista2model _ref2015mod 4 \"--energy ref2015\"\n";
	print "\t2-	take the stuctures relaxed, make a new list mini8newlist and generate loos using rcd\n";
	print "\t3-	time refine_loop.pl ../mini8newlist.txt modeling _ref2015mod 4 \"--energy ref2015\"\n";
	print "\nOLD VERSION MODELING:\n";
	print "\ttime refine_loop.pl ../mini8.txt modeling_old _ref2015mod 4 \"--energy ref2015\"\n";
	print "\nHELP: *Loops can be generated using, for example:\n\ttime mpirun.openmpi -np 4 rcd_mpi_gnu mini8.txt -n 10000 -t 0.99 -r --linear -d 0.5 -x dunbrack.bin -e 5 --energy_file ~/Korp6Dv1/korp6Dv1.bin --bench --loco_best 100 -o mini8_10k1H\n";
	print "\ttime mpirun.openmpi -np 4 rcd_mpi_gnu seok12.txt -n 10000 -t 0.9 -x ~/RCDprot_v103/dunbrack.bin -r --linear -d 0.5 -e 5 --energy_file ~/Korp6Dv1/korp6Dv1.bin --loco_best 100 --bench -o seok12_t90n10k1H\n\n";
	exit;
}
# Load the selected scenario and optional PyRosetta arguments.
my $debug = 0;
my $prog = "refine_loop";
my $rcdfile = $ARGV[0];
my $scenario = $ARGV[1];
my $outsuff = $ARGV[2];
my $np = $ARGV[3];
$optmin .= $ARGV[4]; # Concatenate extra options
$optmin2 .= $ARGV[4]; # Concatenate extra options

# Read the RCD loop list and process each listed loop.
my $data = parseTextFile($rcdfile);

my $nloop = 0;
my $base; # required outside
foreach my $loop (@{$data}) # Screen Loops
{
	print "\n";
	printf "$prog> $scenario SCENARIO. PROCESSING loop %3d: @{$loop}\n",$nloop+1;
	print "\n";
	my $pdb = $loop->[0];
	my $start = $loop->[1];
	my $end = $loop->[2];
	my $chain = $loop->[3];
	my $start2 = $start - 1;
	my $end2 = $end + 1;
	$base = basename($pdb,".pdb");

	
	if($scenario eq "native")
	{
		my $cmd = "time mpirun -np $np $minimizer --base_in $base --ri $start2 --rf $end2 --chain $chain --base_out $base$outsuff --pdb_in $inpref$pdb $optmin";

                system($cmd);
	}
	elsif($scenario eq "modeling")
	{
		
		# Insert RCD's first loop model because Rosetta requires a closed loop in the input.
		my $cmd = "$insert $base"."_loop1.pdb $inpref$pdb $base$outsuff"."_min.pdb noanchors";
		`$cmd` if $run;

		# Refine the structure; the configured neighborhood threshold controls repacking.
		$cmd = "time mpirun -n $np $minimizer --base_in $base --ri $start2 --rf $end2 --chain $chain --pdb_in $base$outsuff"."_min.pdb --base_out $base$outsuff $optmin2";

                system($cmd);
	}
	elsif($scenario eq "modeling_old")
	{
		# Relax the input structure without its native loop before reinserting an RCD model.
		my $cmd = "time $relax --base_in $base --ri $start2 --rf $end2 --chain $chain --pdb_in $inpref$pdb --base_out $base$outsuff"."_noloop $optrelax";

                system($cmd);
		
		# Insert RCD's first loop model because Rosetta requires a closed loop in the input.
		$cmd = "$insert $base"."_loop1.pdb $base$outsuff"."_noloop_min.pdb $base$outsuff"."_min.pdb noanchors";
		`$cmd` if $run;

		# Refine the structure; the configured neighborhood threshold controls repacking.
		$cmd = "time mpirun -n $np $minimizer --base_in $base --ri $start2 --rf $end2 --chain $chain --pdb_in $base$outsuff"."_min.pdb --base_out $base$outsuff $optmin2";

                system($cmd);
	}
	elsif($scenario eq "crista2model")
	{
		# Relax the input structure without its native loop before reinserting an RCD model.
		my $cmd = "time $relax --base_in $base --ri $start2 --rf $end2 --chain $chain --pdb_in $inpref$pdb --base_out $base$outsuff"."_noloop $optrelax";

                system($cmd);
	}
	elsif($scenario eq "crista2model_wl")
	{
		# Relax the input structure while retaining the modeled loop.
		my $cmd = "time $relax --base_in $base --ri $start2 --rf $end2 --chain $chain --pdb_in $inpref$pdb --base_out $base$outsuff $optrelax_withloop";

                system($cmd);
	}
	elsif($scenario eq "stats")
	{
		
	}
	else
	{
		print "$prog> Error! Unknown scenario: $scenario  (Forcing exit!)\n";
		exit;
	}



	unless ($scenario eq "crista2model" || $scenario eq "crista2model_wl"){
		# Combine per-process score and structure outputs, then sort the score table.
		my $kk = "mv $base$outsuff"."_1_rosetta.txt $base$outsuff"."_rosetta.txt";
		`$kk` if $run;
		$kk = "mv $base$outsuff"."_1_loops.pdb"." $base$outsuff"."_loops.pdb";
		`$kk` if $run;
		for(my $i=2; $i<=$np; $i++)
		{
			$kk = "grep -v \"^#\" $base$outsuff"."_$i"."_rosetta.txt >> $base$outsuff"."_rosetta.txt; cat $base$outsuff"."_$i"."_loops.pdb >> $base$outsuff"."_loops.pdb; rm $base$outsuff"."_$i"."_rosetta.txt; rm $base$outsuff"."_$i"."_loops.pdb";
			`$kk` if $run;
		}
		$kk = "head -n 2 $base$outsuff"."_rosetta.txt > $base$outsuff"."_rosettaS.txt; grep -v ^# $base$outsuff"."_rosetta.txt | sort -n -k $col_rosettaE >> $base$outsuff"."_rosettaS.txt;";
		`$kk` if $run;
	
		# Summarize the five highest-ranked candidates.
		$kk = "head -7 $base$outsuff"."_rosettaS.txt |tail -5 >temp; statext.pl temp 1 18 1 all 2 > $base$outsuff"."_aveTop5.txt;rm temp ";
		`$kk` if $run;
		
		
		print "\n";
		print "Refined PDB is:  $base$outsuff"."_loops.pdb\n";
		print "Refined Rosetta: $base$outsuff"."_rosetta.txt\n";
		print "Sorted Rosetta:  $base$outsuff"."_rosettaS.txt\n";
		print "Average Top5:  $base$outsuff"."_aveTop5.txt\n";
		print "\n";
	}
}



unless ($scenario eq "crista2model" || $scenario eq "crista2model_wl"){
	# Summarize the average score across each structure’s five top candidates.
	my $kk = "head -n 1 $base$outsuff"."_rosetta.txt > rosettaAve5$outsuff.txt; for i in *$outsuff"."_aveTop5.txt; do var=`tail -3 \$i|head -1` ; echo -e \"\${i:0:20} \${var:8}\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosettaAve5$outsuff.txt";
	open my $fh, '-|', qw(/bin/bash -c), $kk;
	close($fh);
	sleep(1);
	# Summarize each structure’s best candidate.
	$kk = "head -n 1 $base$outsuff"."_rosetta.txt > rosetta1$outsuff.txt; for i in *$outsuff"."_rosettaS.txt; do var=\`awk \'NR==3\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosetta1$outsuff.txt";
	open $fh, '-|', qw(/bin/bash -c), $kk;
	close($fh);
	sleep(1);
	# Summarize the best candidate among each structure’s top two.
	$kk = "head -n 1 $base$outsuff"."_rosetta.txt > rosetta2$outsuff.txt; for i in *$outsuff"."_rosettaS.txt; do var=\`awk \'NR>=3 && NR<3+2\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosetta2$outsuff.txt";
	open $fh, '-|', qw(/bin/bash -c), $kk;
	close($fh);
	sleep(1);
	# Summarize the best candidate among each structure’s top five.
	$kk = "head -n 1 $base$outsuff"."_rosetta.txt > rosetta5$outsuff.txt; for i in *$outsuff"."_rosettaS.txt; do var=\`awk \'NR>=3 && NR<3+5\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosetta5$outsuff.txt";
	open $fh, '-|', qw(/bin/bash -c), $kk;
	close($fh);
	sleep(1);
	# Summarize the best candidate among each structure’s top ten.
	$kk = "head -n 1 $base$outsuff"."_rosetta.txt > rosetta10$outsuff.txt; for i in *$outsuff"."_rosettaS.txt; do var=\`awk \'NR>=3 && NR<3+10\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosetta10$outsuff.txt";
	open $fh, '-|', qw(/bin/bash -c), $kk;
	close($fh);

	sleep(1);
	# Summarize the best candidate among each structure’s top twenty.
	$kk = "head -n 1 $base$outsuff"."_rosetta.txt > rosetta20$outsuff.txt; for i in *$outsuff"."_rosettaS.txt; do var=\`awk \'NR>=3 && NR<3+20\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosetta20$outsuff.txt; rm temp";
	open $fh, '-|', qw(/bin/bash -c), $kk;
	close($fh);
	sleep(1);
	print "\n$prog> That's all folks!\n";
	# Summarize the best candidate among each structure’s top fifty.
	$kk = "head -n 1 $base$outsuff"."_rosetta.txt > rosetta50$outsuff.txt; for i in *$outsuff"."_rosettaS.txt; do var=\`awk \'NR>=3 && NR<3+50\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosetta50$outsuff.txt; rm temp";
	open $fh, '-|', qw(/bin/bash -c), $kk;
	close($fh);
	sleep(1);
	# Summarize the best candidate among each structure’s top one hundred.
	$kk = "head -n 1 $base$outsuff"."_rosetta.txt > rosetta100$outsuff.txt; for i in *$outsuff"."_rosettaS.txt; do var=\`awk \'NR>=3 && NR<3+100\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosetta100$outsuff.txt; rm temp";
	open $fh, '-|', qw(/bin/bash -c), $kk;
	close($fh);
	sleep(1);
	exit;
}


# Parse each non-comment row of the RCD loop list into its fields.
sub parseTextFile
{
	my $file = shift; # Input file name
	open(IN,$file) or die "\nFailed to open $file\n";

	my @data = ();
	my $nline = 0; # Line index (counter)
	print "Reading loop list: $file\n";
	while(<IN>)
	{
		next if /^#/; # reading only non-# begining lines
		my @line = split(/\s+/, $_);
		push(@data,\@line);
		$nline++;
	}
	print "Found $nline loop(s).\n";
	close IN;
	
	return \@data; # Return reference to data
}


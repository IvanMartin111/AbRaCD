#!/usr/bin/perl
#
# RCD+ loop refinement for "native" or "modeling" scenarios
#
# by Mon (05/02/2019)
#

use strict;
use Storable;
#use Expect;
use File::Basename;



# SOME DEFINITIONS

# Path to the Rosetta's loop minimization python (e.g. ../loops_energyXX.py). Relative to the execution directory.
my $minimizer = "loops_energy13.py";
my $relax = "monrelax4.py";
my $insert = "insertpdb.pl";
my $inpref = "../";
my $col_rosettaE = 9; # Column index (for -k option of sort command)
my $col_rmsd = 12; # Column index of PyRosetta refined RMSD
my $run = 1; # =0, Disables execution and forces plain text dump of commands; otherwise run.

# PyRosetta options for minimization
# NATIVE scenario:
my $optmin = "--nmin 5 --initneighbor -1 --minneighbor 0 --neighbor -1 --irepack 0 --mintype dfpmin_armijo_nonmonotone --nrep 1 --fas=0.05,0.30,0.60,1.00 --ntpsi --verbose error ";
# MODELING scenario:
my $optrelax = "--nmin 5 --initneighbor 10 --nminneigh 100 --irepack 0 --mintype dfpmin_armijo_nonmonotone --fas=0.05,0.30,0.60,1.00 --ntpsi --deleteloop ";
my $optrelax_withloop = "--nmin 5 --initneighbor 10 --nminneigh 100 --irepack 0 --mintype dfpmin_armijo_nonmonotone --ntpsi";
my $optmin2 = "--nmin 5 --initneighbor -1 --minneighbor 5 --neighbor -1 --irepack 0 --mintype dfpmin_armijo_nonmonotone --nrep 1 --fas=0.05,0.30,0.60,1.00 --ntpsi --verbose error ";


# INPUT PARSER
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
	print "\ttime refine_loops.pl ../mini8.txt native _ref2015nat 4 \"--energy ref2015\"\n";
	print "\ttime refine_loops.pl ../mini8.txt modeling _ref2015mod 4 \"--energy ref2015\"\n";
	print "\ttime refine_loops.pl ../mini8.txt stats _ref2015mod 4 \"--energy ref2015\"\t#scenario for only recalculate thet stats files without using rosetta\n";
	print "\nCRISTALTOMODEL steps:\n";
	print "\t1-	time refine_loops.pl ../mini8.txt crista2model _ref2015mod 4 \"--energy ref2015\"\n";
	print "\t2-	take the stuctures relaxed, make a new list mini8newlist and generate loos using rcd\n";
	print "\t3-	time refine_loops.pl ../mini8newlist.txt modeling _ref2015mod 4 \"--energy ref2015\"\n";
	print "\nOLD VERSION MODELING:\n";
	print "\ttime refine_loops.pl ../mini8.txt modeling_old _ref2015mod 4 \"--energy ref2015\"\n";
	print "\nHELP: *Loops can be generated using, for example:\n\ttime mpirun.openmpi -np 4 rcd_mpi_gnu mini8.txt -n 10000 -t 0.99 -r --linear -d 0.5 -x dunbrack.bin -e 5 --energy_file ~/Korp6Dv1/korp6Dv1.bin --bench --loco_best 100 -o mini8_10k1H\n";
	print "\ttime mpirun.openmpi -np 4 rcd_mpi_gnu seok12.txt -n 10000 -t 0.9 -x ~/RCDprot_v103/dunbrack.bin -r --linear -d 0.5 -e 5 --energy_file ~/Korp6Dv1/korp6Dv1.bin --loco_best 100 --bench -o seok12_t90n10k1H\n\n";
	exit;
}
# Some commands:
# for i in *_ref2015nat_rosettaS.txt; do kk=`awk 'NR==3' $i`; echo -e "${i:0:15} $kk"; done > rosetta_ref2015nat.txt
# for i in *_ref2015mod_rosettaS.txt; do kk=`awk 'NR==3' $i`; echo -e "${i:0:15} $kk"; done > rosetta_ref2015mod.txt

# Loading things        
my $debug = 0;
my $prog = "refine_loops";
my $rcdfile = $ARGV[0];
my $scenario = $ARGV[1];
my $outsuff = $ARGV[2];
my $np = $ARGV[3];
$optmin .= $ARGV[4]; # Concatenate extra options
$optmin2 .= $ARGV[4]; # Concatenate extra options

# Read input file
my $data = parseTextFile($rcdfile);

my $nloop = 0;
my $base; # required outside
foreach my $loop (@{$data}) # Screen Loops
{
	printf "$prog> $scenario SCENARIO. PROCESSING loop %3d: @{$loop}\n",$nloop+1;
	my $pdb = $loop->[0];
	my $start = $loop->[1];
	my $end = $loop->[2];
	my $chain = $loop->[3];
	my $start2 = $start - 1;
	my $end2 = $end + 1;
	$base = basename($pdb,".pdb");

	
	# MAIN LOOP for NATIVE scenario
	if($scenario eq "native")
	{
		my $cmd = "time mpirun -np $np $minimizer --base_in $base --ri $start2 --rf $end2 --chain $chain --base_out $base$outsuff --pdb_in $inpref$pdb $optmin";

		# Using "Expect.pm" to dump program output during command execution
		# Install Expect.pm using: perl -MCPAN -e 'install Bundle::Expect'
		print "EXECUTING COMMAND: $cmd\n";
                system($cmd);
#		if($run)
#		{
#			my $process = Expect->spawn($cmd) or die "Cannot spawn $cmd: $!\n";
#			while ($process->expect(undef))
#			{
#			  print $process->before();
#			}
#		}
	}
	# MAIN LOOP for MODELING scenario
	elsif($scenario eq "modeling")
	{
#		# 5. Main PDB structure must be relaxed first in Rosetta potential without the native loop (--deleteloop).
#		#    (This stage can be carried out in parallel with the next one.)
#		my $cmd = "time $relax --base_in $base --ri $start2 --rf $end2 --chain $chain --pdb_in $inpref$pdb --base_out $base$outsuff"."_noloop $optrelax";
#
#		# Using "Expect.pm" to dump program output during command execution
#		# Install Expect.pm using: perl -MCPAN -e 'install Bundle::Expect'
#		print "EXECUTING COMMAND: $cmd\n";
#                system($cmd);
#		if($run)
#		{
#			my $process = Expect->spawn($cmd) or die "Cannot spawn $cmd: $!\n";
#			while ($process->expect(undef))
#			{
#			  print $process->before();
#			}
#		}
		
		# 7. Rosetta requires that some closed loop be present in the complete input structure.
		#    (native loop, if present, was deleted by --deleteloop option).
		#    The script "insertpdb.pl" will insert the first closed loop, "*_loop1.pdb",
		#    that was automatically generated by RCD into the minimized input PDB.
		my $cmd = "$insert $base"."_loop1.pdb $inpref$pdb $base$outsuff"."_min.pdb noanchors";
		print "EXECUTING COMMAND: $cmd\n";
		`$cmd` if $run;

		# 8. Run Rosetta refinement in PyRosetta. A 5A threshold will trigger repacking.
		#    Use the "input_withloop1.pdb" instead of "input.pdb" if neccessary.
		$cmd = "time mpirun -n $np $minimizer --base_in $base --ri $start2 --rf $end2 --chain $chain --pdb_in $base$outsuff"."_min.pdb --base_out $base$outsuff $optmin2";

		# Using "Expect.pm" to dump program output during command execution
		# Install Expect.pm using: perl -MCPAN -e 'install Bundle::Expect'
		print "EXECUTING COMMAND: $cmd\n";
                system($cmd);
#		if($run)
#		{
#			my $process = Expect->spawn($cmd) or die "Cannot spawn $cmd: $!\n";
#			while ($process->expect(undef))
#			{
#			  print $process->before();
#			}
#		}
	}
	# MAIN LOOP for MODELING scenario
	elsif($scenario eq "modeling_old")
	{
		# 5. Main PDB structure must be relaxed first in Rosetta potential without the native loop (--deleteloop).
		#    (This stage can be carried out in parallel with the next one.)
		my $cmd = "time $relax --base_in $base --ri $start2 --rf $end2 --chain $chain --pdb_in $inpref$pdb --base_out $base$outsuff"."_noloop $optrelax";

		# Using "Expect.pm" to dump program output during command execution
		# Install Expect.pm using: perl -MCPAN -e 'install Bundle::Expect'
		print "EXECUTING COMMAND: $cmd\n";
                system($cmd);
#		if($run)
#		{
#			my $process = Expect->spawn($cmd) or die "Cannot spawn $cmd: $!\n";
#			while ($process->expect(undef))
#			{
#			  print $process->before();
#			}
#		}
		
		# 7. Rosetta requires that some closed loop be present in the complete input structure.
		#    (native loop, if present, was deleted by --deleteloop option).
		#    The script "insertpdb.pl" will insert the first closed loop, "*_loop1.pdb",
		#    that was automatically generated by RCD into the minimized input PDB.
		$cmd = "$insert $base"."_loop1.pdb $base$outsuff"."_noloop_min.pdb $base$outsuff"."_min.pdb noanchors";
		print "EXECUTING COMMAND: $cmd\n";
		`$cmd` if $run;

		# 8. Run Rosetta refinement in PyRosetta. A 5A threshold will trigger repacking.
		#    Use the "input_withloop1.pdb" instead of "input.pdb" if neccessary.
		$cmd = "time mpirun -n $np $minimizer --base_in $base --ri $start2 --rf $end2 --chain $chain --pdb_in $base$outsuff"."_min.pdb --base_out $base$outsuff $optmin2";

		# Using "Expect.pm" to dump program output during command execution
		# Install Expect.pm using: perl -MCPAN -e 'install Bundle::Expect'
		print "EXECUTING COMMAND: $cmd\n";
                system($cmd);
#		if($run)
#		{
#			my $process = Expect->spawn($cmd) or die "Cannot spawn $cmd: $!\n";
#			while ($process->expect(undef))
#			{
#			  print $process->before();
#			}
#		}
	}
	elsif($scenario eq "crista2model")
	{
		# 5. Main PDB structure must be relaxed first in Rosetta potential without the native loop (--deleteloop).
		#    (This stage can be carried out in parallel with the next one.)
		my $cmd = "time $relax --base_in $base --ri $start2 --rf $end2 --chain $chain --pdb_in $inpref$pdb --base_out $base$outsuff"."_noloop $optrelax";

		# Using "Expect.pm" to dump program output during command execution
		# Install Expect.pm using: perl -MCPAN -e 'install Bundle::Expect'
		print "EXECUTING COMMAND: $cmd\n";
                system($cmd);
#		if($run)
#		{
#			my $process = Expect->spawn($cmd) or die "Cannot spawn $cmd: $!\n";
#			while ($process->expect(undef))
#			{
#			  print $process->before();
#			}
#		}
	}
	elsif($scenario eq "crista2model_wl")
	{
		# 5. Main PDB structure must be relaxed first in Rosetta potential without the native loop (--deleteloop).
		#    (This stage can be carried out in parallel with the next one.)
		my $cmd = "time $relax --base_in $base --ri $start2 --rf $end2 --chain $chain --pdb_in $inpref$pdb --base_out $base$outsuff $optrelax_withloop";

		# Using "Expect.pm" to dump program output during command execution
		# Install Expect.pm using: perl -MCPAN -e 'install Bundle::Expect'
		print "EXECUTING COMMAND: $cmd\n";
                system($cmd);
#		if($run)
#		{
#			my $process = Expect->spawn($cmd) or die "Cannot spawn $cmd: $!\n";
#			while ($process->expect(undef))
#			{
#			  print $process->before();
#			}
#		}
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
		# Merge all *_?_rosetta.txt and *_?_loops.pdb files
		# Merge _rosetta.txt files
		my $kk = "mv $base$outsuff"."_1_rosetta.txt $base$outsuff"."_rosetta.txt";
		print "EXECUTING: $kk\n";
		`$kk` if $run;
		$kk = "mv $base$outsuff"."_1_loops.pdb"." $base$outsuff"."_loops.pdb";
		print "EXECUTING: $kk\n";
		`$kk` if $run;
		for(my $i=2; $i<=$np; $i++)
		{
			$kk = "grep -v \"^#\" $base$outsuff"."_$i"."_rosetta.txt >> $base$outsuff"."_rosetta.txt; cat $base$outsuff"."_$i"."_loops.pdb >> $base$outsuff"."_loops.pdb; rm $base$outsuff"."_$i"."_rosetta.txt; rm $base$outsuff"."_$i"."_loops.pdb";
			print "kk= $kk\n";
			`$kk` if $run;
		}
		$kk = "head -n 2 $base$outsuff"."_rosetta.txt > $base$outsuff"."_rosettaS.txt; grep -v ^# $base$outsuff"."_rosetta.txt | sort -n -k $col_rosettaE >> $base$outsuff"."_rosettaS.txt;";
		print "EXECUTING: $kk\n";
		`$kk` if $run;
	
		print "Refined PDB is:  $base$outsuff"."_loops.pdb\n";
		print "Refined Rosetta: $base$outsuff"."_rosetta.txt\n";
		print "Sorted Rosetta:  $base$outsuff"."_rosettaS.txt\n";
	
		#statext for mean of top5
		$kk = "head -7 $base$outsuff"."_rosettaS.txt |tail -5 >temp; statext.pl temp 1 18 1 all 2 > $base$outsuff"."_aveTop5.txt;rm temp ";
		print "EXECUTING: $kk\n";
		`$kk` if $run;
		
		
		print "Average Top5:  $base$outsuff"."_aveTop5.txt\n";
	}
}



unless ($scenario eq "crista2model" || $scenario eq "crista2model_wl"){
	# Get Average-5 energy decoys statistics
	my $kk = "head -n 1 $base$outsuff"."_rosetta.txt > rosettaAve5$outsuff.txt; for i in *$outsuff"."_aveTop5.txt; do var=`tail -3 \$i|head -1` ; echo -e \"\${i:0:20} \${var:8}\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosettaAve5$outsuff.txt";
	print "EXECUTING: $kk\n";
	open my $fh, '-|', qw(/bin/bash -c), $kk;
	close($fh);
	sleep(1);
	# Get best-1 energy decoys statistics
	$kk = "head -n 1 $base$outsuff"."_rosetta.txt > rosetta1$outsuff.txt; for i in *$outsuff"."_rosettaS.txt; do var=\`awk \'NR==3\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosetta1$outsuff.txt";
	print "EXECUTING: $kk\n";
	open $fh, '-|', qw(/bin/bash -c), $kk;
	close($fh);
	sleep(1);
	# Get best-2 energy decoys statistics
	$kk = "head -n 1 $base$outsuff"."_rosetta.txt > rosetta2$outsuff.txt; for i in *$outsuff"."_rosettaS.txt; do var=\`awk \'NR>=3 && NR<3+2\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosetta2$outsuff.txt";
	print "EXECUTING: $kk\n";
	open $fh, '-|', qw(/bin/bash -c), $kk;
	close($fh);
	sleep(1);
	# Get best-5 energy decoys statistics
	$kk = "head -n 1 $base$outsuff"."_rosetta.txt > rosetta5$outsuff.txt; for i in *$outsuff"."_rosettaS.txt; do var=\`awk \'NR>=3 && NR<3+5\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosetta5$outsuff.txt";
	print "EXECUTING: $kk\n";
	open $fh, '-|', qw(/bin/bash -c), $kk;
	close($fh);
	sleep(1);
	# Get best-10 energy decoys statistics
	$kk = "head -n 1 $base$outsuff"."_rosetta.txt > rosetta10$outsuff.txt; for i in *$outsuff"."_rosettaS.txt; do var=\`awk \'NR>=3 && NR<3+10\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosetta10$outsuff.txt";
	print "EXECUTING: $kk\n";
	open $fh, '-|', qw(/bin/bash -c), $kk;
	close($fh);

	sleep(1);
	# Get best-20 energy decoys statistics
	$kk = "head -n 1 $base$outsuff"."_rosetta.txt > rosetta20$outsuff.txt; for i in *$outsuff"."_rosettaS.txt; do var=\`awk \'NR>=3 && NR<3+20\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosetta20$outsuff.txt; rm temp";
	print "EXECUTING: $kk\n";
	open $fh, '-|', qw(/bin/bash -c), $kk;
	close($fh);
	sleep(1);
	print "$prog> That's all folks!\n";
	# Get best-50 energy decoys statistics
	$kk = "head -n 1 $base$outsuff"."_rosetta.txt > rosetta50$outsuff.txt; for i in *$outsuff"."_rosettaS.txt; do var=\`awk \'NR>=3 && NR<3+50\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosetta50$outsuff.txt; rm temp";
	print "EXECUTING: $kk\n";
	open $fh, '-|', qw(/bin/bash -c), $kk;
	close($fh);
	sleep(1);
	# Get best-100 energy decoys statistics
	$kk = "head -n 1 $base$outsuff"."_rosetta.txt > rosetta100$outsuff.txt; for i in *$outsuff"."_rosettaS.txt; do var=\`awk \'NR>=3 && NR<3+100\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosetta100$outsuff.txt; rm temp";
	print "EXECUTING: $kk\n";
	open $fh, '-|', qw(/bin/bash -c), $kk;
	close($fh);
	sleep(1);
	exit;
}


#################################################################
# FUNCTIONS
#################################################################

# Read text file line by line and separate each line into fields
sub parseTextFile
{
	my $file = shift; # Input file name
	open(IN,$file) or die "\nFailed to open $file\n";

	my @data = ();
	my $nline = 0; # Line index (counter)
	print "parseTextFile> Reading file: $file --> ";
	while(<IN>)
	{
		next if /^#/; # reading only non-# begining lines
		my @line = split(/\s+/, $_);
		# print "@line\n";
		push(@data,\@line);
		$nline++;
	}
	print " $nline lines\n";
	close IN;
	
	return \@data; # Return reference to data
}


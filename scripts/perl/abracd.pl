#!/usr/bin/env perl

# Provisional end-to-end launcher for AbRaCD.
use strict;
use warnings;
use FindBin qw($Bin);
use File::Basename qw(basename);
use File::Spec;

my $PREPROCESS = File::Spec->catfile($Bin, 'AbRaCD_chothia_preprocces.pl');
my $REFINE     = File::Spec->catfile($Bin, 'refine_loops.pl');
my $REPO_ROOT  = File::Spec->rel2abs(File::Spec->catdir($Bin, '..', '..'));
my $MAPS_DIR   = File::Spec->catdir($REPO_ROOT, 'data', 'Maps');
my $CDR_MAP    = File::Spec->catfile($MAPS_DIR, 'CDRS.bin');
my $KORP_MAP   = File::Spec->catfile($MAPS_DIR, 'korp6Dv1.bin');
my $DUN_MAP    = File::Spec->catfile($MAPS_DIR, 'dunbrack.bin');

# Default RCD settings, based on scripts/bash/script_allshort_CDRs_RCD.sh.
my $NP       = 1;
my %CDR_ARGS = (
    H1 => ['--H1', 0.5, 0.2, 20000, 2000], H2 => ['--H2', 0.9, 0.5, 20000, 2000],
    H3 => ['--H3', 0.5, 0.5, 50000, 1000], L1 => ['--L1', 0.9, 0.5, 20000, 2000],
    L2 => ['--L2', 0.9, 0.5, 20000, 2000], L3 => ['--L3', 0.9, 0.5, 20000, 2000],
);
my %CDR_EXTRA_ARGS = (
    H3 => ['--kink_merge', 0.9, '--kink_mergeA', 0.9, '--tkink', 0.70, '--tkinkA', 0.70],
);

sub usage {
    my ($fh) = @_;
    print {$fh} <<'USAGE';
abracd.pl — run the AbRaCD workflow for one CDR

Usage:
  abracd.pl [STRUCTURE.pdb | LIST.txt] [--raw] [--purge] --CDR CDR [--NP X] [--NLOOPS X] [--LOCOBEST X] [--repRCD X] [--resROS Y]
  abracd.pl [--raw] --CDR CDR [--NP X] [--NLOOPS X] [--LOCOBEST X] [--repRCD X] [--resROS Y]

Input:
  No structure/list:   process all PDB files in the current directory.
  STRUCTURE.pdb:       process one structure.
  LIST.txt:            one PDB per line (relative paths use the list directory).
  --raw:               skip Chothia preprocessing; place immediately after
                       the input argument (or first when processing all PDBs).
  --purge:             use the same input file as the purge file. It must have
                       one row per structure: PDB  H_chain  L_chain (use - for
                       a missing chain). Only applies without --raw.
  --CDR CDR:           required. Accepted values: H1, H2, H3, L1, L2, L3.
  --NP X:              use the MPI RCD executable with X processes (default: serial RCD).
  --NLOOPS X:          RCD loop sampling count (default: 50000 for H3; 20000 otherwise).
  --LOCOBEST X:        RCD loco_best count (default: 1000 for H3; 2000 otherwise).
  --repRCD X:          number of RCD runs (default: 1).
  --resROS Y:          Rosetta refinements per RCD run (default: 1).

Examples:
  abracd.pl antibody.pdb --CDR H3
  abracd.pl structures.txt --CDR H3 --repRCD 3 --resROS 2
  abracd.pl purge.txt --purge --CDR L3
  abracd.pl structures.txt --raw --CDR H3
  abracd.pl --raw --CDR H3

If a step fails, the failed command and an example invocation are shown.
USAGE
}

sub fail {
    my ($message) = @_;
    print STDERR "ERROR: $message\n\n";
    usage(*STDERR);
    exit 2;
}

sub run_command {
    my (@cmd) = @_;
    print "\n> ", join(' ', map { /\s/ ? "'$_'" : $_ } @cmd), "\n";
    my $status = system { $cmd[0] } @cmd;
    return if $status == 0;
    my $detail = $status == -1 ? "could not execute: $!"
                : ($status & 127) ? "terminated by signal " . ($status & 127)
                : "exit code " . ($status >> 8);
    fail("Step failed ($detail). Command: " . join(' ', @cmd)
       . "\nExample: abracd.pl [STRUCTURE.pdb|LIST.txt] --CDR H3");
}

sub find_executable {
    my ($name) = @_;
    return $name if File::Spec->file_name_is_absolute($name) && -x $name;
    for my $dir (File::Spec->path()) {
        my $candidate = File::Spec->catfile($dir, $name);
        return $candidate if -f $candidate && -x $candidate;
    }
    return;
}

sub check_python_modules {
    my ($python, @modules) = @_;
    my $code = 'import importlib.util, sys; missing=[m for m in sys.argv[1:] if importlib.util.find_spec(m) is None]; print(" ".join(missing)); sys.exit(bool(missing))';
    my $status = system { $python } $python, '-c', $code, @modules;
    return $status == 0;
}

my @args = @ARGV;
if (@args == 1 && ($args[0] eq '-h' || $args[0] eq '--help')) { usage(*STDOUT); exit 0; }

my $input;
if (@args && $args[0] !~ /^--/) { $input = shift @args; }
my $raw = 0;
if (@args && $args[0] eq '--raw') { $raw = 1; shift @args; }
elsif (grep { $_ eq '--raw' } @args) {
    fail('--raw must appear immediately after the structure/list argument (or first when processing all PDB files).');
}

my $purge = 0;
my ($cdr, $rep_rcd, $res_ros, $nloops_override, $locobest_override) = (undef, 1, 1, undef, undef);
my $np_requested = 0;
while (@args) {
    my $arg = shift @args;
    if ($arg eq '--purge') { $purge = 1; }
    elsif ($arg eq '--CDR') { $cdr = shift @args; }
    elsif ($arg eq '--NP') { $NP = shift @args; $np_requested = 1; }
    elsif ($arg eq '--NLOOPS') { $nloops_override = shift @args; }
    elsif ($arg eq '--LOCOBEST') { $locobest_override = shift @args; }
    elsif ($arg eq '--repRCD') { $rep_rcd = shift @args; }
    elsif ($arg eq '--resROS') { $res_ros = shift @args; }
    else { fail("Unknown argument or invalid order: $arg"); }
}
fail('--CDR is required. For example, specify --CDR H3.') unless defined $cdr;
$cdr = uc $cdr;
fail("Invalid CDR: $cdr. Accepted values: H1, H2, H3, L1, L2, L3.") unless exists $CDR_ARGS{$cdr};
fail('--purge and --raw cannot be used together.') if $purge && $raw;
for my $pair (['--repRCD', $rep_rcd], ['--resROS', $res_ros]) {
    fail("$pair->[0] must be a positive integer.")
        unless defined($pair->[1]) && $pair->[1] =~ /^\d+$/ && $pair->[1] > 0;
}
for my $pair (['--NP', $NP], ['--NLOOPS', $nloops_override], ['--LOCOBEST', $locobest_override]) {
    next unless defined $pair->[1];
    fail("$pair->[0] must be a positive integer.")
        unless $pair->[1] =~ /^\d+$/ && $pair->[1] > 0;
}
fail("Preprocessor not found: $PREPROCESS") unless -x $PREPROCESS || -f $PREPROCESS;
fail("Refinement script not found: $REFINE") unless -f $REFINE;
for my $map ($CDR_MAP, $KORP_MAP, $DUN_MAP) {
    fail("Required RCD map not found: $map") unless -f $map;
}
my $python = find_executable('python3') || find_executable('python');
fail('Python 3 was not found on PATH; PyRosetta and mpi4py cannot be checked.') unless $python;
my @missing_modules;
push @missing_modules, 'PyRosetta (python module: pyrosetta)' unless check_python_modules($python, 'pyrosetta');
push @missing_modules, 'mpi4py (python module: mpi4py)' unless check_python_modules($python, 'mpi4py');
fail("Required Python dependencies are missing for '$python':\n  - " . join("\n  - ", @missing_modules)) if @missing_modules;

my $rcd_executable;
my $use_mpi = 0;
if ($np_requested) {
    $rcd_executable = find_executable('rcd_mpi_gcc')
        || find_executable('rcd_mpi_intel') || find_executable('rcd_mpi');
    fail('--NP was requested, but no MPI RCD executable was found on PATH (looked for rcd_mpi_gcc, rcd_mpi_intel, rcd_mpi).')
        unless $rcd_executable;
    $use_mpi = 1;
    my $mpirun = find_executable('mpirun');
    fail('--NP was requested, but mpirun was not found on PATH.') unless $mpirun;
} else {
    $rcd_executable = find_executable('rcd_gcc') || find_executable('rcd');
    unless ($rcd_executable) {
        $rcd_executable = find_executable('rcd_mpi_gcc')
            || find_executable('rcd_mpi_intel') || find_executable('rcd_mpi');
        fail('No RCD executable was found on PATH (looked for rcd_gcc, rcd, rcd_mpi_gcc, rcd_mpi_intel, rcd_mpi).')
            unless $rcd_executable;
        $use_mpi = 1;
        my $mpirun = find_executable('mpirun');
        fail('Only an MPI RCD executable is available, but mpirun was not found on PATH.') unless $mpirun;
    }
}
my ($flag, $nokink, $dval, $default_nloops, $default_locobest) = @{$CDR_ARGS{$cdr}};
my $nloops = defined($nloops_override) ? $nloops_override : $default_nloops;
my $locobest = defined($locobest_override) ? $locobest_override : $default_locobest;
my @cdr_extra_args = @{$CDR_EXTRA_ARGS{$cdr} || []};
print "Preflight passed. Maps found; Python dependencies available; RCD executable: $rcd_executable\n";
print $use_mpi ? "RCD mode: MPI ($NP processes)\n" : "RCD mode: serial (default NP=1)\n";
fail('With --purge, provide a .txt list with rows in this format: PDB H_chain L_chain.')
    if $purge && (!defined($input) || $input !~ /\.txt$/i);

my @pre_cmd = ($^X, $PREPROCESS);
push @pre_cmd, '--purge', $input if $purge;
push @pre_cmd, $input if defined($input) && !$purge;
run_command(@pre_cmd) unless $raw;

my $rcd_input = "AB_loop_AbRaCD_${cdr}.txt";
fail("RCD input file '$rcd_input' not found. In --raw mode, create this file before running; it must contain the RCD entries for CDR $cdr.")
    unless -f $rcd_input;
for my $rep (1 .. $rep_rcd) {
    my $out = "abracd_${cdr}_rep${rep}";
    my @rcd_cmd = ($use_mpi ? ('mpirun', '-np', $NP) : (), $rcd_executable, $rcd_input,
        '--ABmap', $CDR_MAP, $flag, @cdr_extra_args, '--nokink_merge', $nokink,
        '--nterm_rama', 0, '-r', '-e', 56, '--energy_file', $KORP_MAP, '--bench',
        '-n', $nloops, '--loco_best', $locobest, '-x', $DUN_MAP,
        '-d', $dval, '--linear', '-t', 0.99, '-o', $out);
    run_command(@rcd_cmd);

    for my $ros (1 .. $res_ros) {
        my $dir = $out;
        fail("RCD completed successfully but did not create directory '$dir'.") unless -d $dir;
        my $old = File::Spec->rel2abs('.');
        chdir $dir or fail("Cannot enter '$dir': $!");
        my $suffix = "_abracd_${cdr}_rep${rep}_ros${ros}";
        my @refine_cmd = ($^X, $REFINE, "../$rcd_input", 'native', $suffix, $NP, '--energy ref2015');
        my $status = system { $refine_cmd[0] } @refine_cmd;
        chdir $old or fail("Cannot return to directory '$old': $!");
        if ($status != 0) {
            fail("Rosetta refinement failed for $out (exit code " . ($status >> 8) . ").\nCommand: " . join(' ', @refine_cmd)
                . "\nExample: abracd.pl " . (defined($input) ? $input : '') . " --CDR $cdr --repRCD $rep_rcd --resROS $res_ros");
        }
    }
}

print "\nWorkflow completed for $cdr. RCD runs: $rep_rcd; Rosetta refinements per RCD run: $res_ros.\n";
exit 0;

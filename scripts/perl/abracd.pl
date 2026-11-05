#!/usr/bin/env perl

# Provisional end-to-end launcher for AbRaCD.
use strict;
use warnings;
use FindBin qw($Bin);
use File::Basename qw(basename);
use File::Spec;

my $PREPROCESS = File::Spec->catfile($Bin, 'AbRaCD_chothia_preprocces.pl');
my $REFINE     = File::Spec->catfile($Bin, 'refine_loops.pl');

# Initial RCD settings, based on scripts/bash/script_allshort_CDRs_RCD.sh.
# These defaults will be made configurable as the workflow is specified.
my $NP       = 10;
my $NLOOPS   = 20000;
my $LOCOBEST = 200;
my $ENERGY   = 'korp3A90c18mr1f10_Allno50_6Dc14mr3r10n36nc8b4nb1yz8bf18_map.bin';
my %CDR_ARGS = (
    H1 => ['--H1', 0.4, 0.2], H2 => ['--H2', 0.9, 0.5],
    H3 => ['--H3', 0.9, 0.5], L1 => ['--L1', 0.9, 0.5],
    L2 => ['--L2', 0.9, 0.5], L3 => ['--L3', 0.9, 0.5],
);

sub usage {
    my ($fh) = @_;
    print {$fh} <<'USAGE';
abracd.pl — run the AbRaCD workflow for one CDR

Usage:
  abracd.pl [STRUCTURE.pdb | LIST.txt] [--raw] [--purge] --CDR CDR [--repRCD X] [--resROS Y]
  abracd.pl [--raw] --CDR CDR [--repRCD X] [--resROS Y]

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
my ($cdr, $rep_rcd, $res_ros) = (undef, 1, 1);
while (@args) {
    my $arg = shift @args;
    if ($arg eq '--purge') { $purge = 1; }
    elsif ($arg eq '--CDR') { $cdr = shift @args; }
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
fail("Preprocessor not found: $PREPROCESS") unless -x $PREPROCESS || -f $PREPROCESS;
fail("Refinement script not found: $REFINE") unless -f $REFINE;
fail('With --purge, provide a .txt list with rows in this format: PDB H_chain L_chain.')
    if $purge && (!defined($input) || $input !~ /\.txt$/i);

my @pre_cmd = ($^X, $PREPROCESS);
push @pre_cmd, '--purge', $input if $purge;
push @pre_cmd, $input if defined($input) && !$purge;
run_command(@pre_cmd) unless $raw;

my $rcd_input = "AB_loop_AbRaCD_${cdr}.txt";
fail("RCD input file '$rcd_input' not found. In --raw mode, create this file before running; it must contain the RCD entries for CDR $cdr.")
    unless -f $rcd_input;
my ($flag, $kink, $dval) = @{$CDR_ARGS{$cdr}};

for my $rep (1 .. $rep_rcd) {
    my $out = "abracd_${cdr}_rep${rep}";
    my @rcd_cmd = ('mpirun', '-np', $NP, 'rcd_mpi', $rcd_input,
        '--ABmap', 'CDRS.bin', $flag, '--nokink_merge', $kink,
        '-r', '-e', 56, '--energy_file', $ENERGY, '--bench',
        '-n', $NLOOPS, '--loco_best', $LOCOBEST, '-x', 'dunbrack.bin',
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

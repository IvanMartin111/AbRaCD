#!/usr/bin/perl
use strict;
use warnings;
use File::Basename qw(dirname basename);
use File::Spec;
use Cwd qw(abs_path);

sub usage {
    print <<'USAGE';
Prepare antibody engineering inputs for RCD

Usage:
  ab_engineering_prep.pl RCD_INPUT.txt

The input must contain standard RCD rows:
  PDB_file  start  end  chain  sequence  length  CDR

Only structures named in the RCD input are processed. For each row, the
specified loop residues are removed from that chain, residues before the loop
keep their numbering, and residues after it are shifted to accommodate the
design sequence length. Other chains and non-coordinate records are copied.

Outputs are written beside the input RCD list and source structures:
  <source>_design.pdb
  AbRaCD_design_<CDR>.txt
USAGE
}

if (@ARGV == 1 && ($ARGV[0] eq '-h' || $ARGV[0] eq '--help')) { usage(); exit 0; }
die "Usage: ab_engineering_prep.pl RCD_INPUT.txt\n" unless @ARGV == 1;
my $input_list = shift @ARGV;
die "RCD input list not found: $input_list\n" unless -f $input_list;
my $list_abs = abs_path($input_list);
my $list_dir = dirname($list_abs);

my %design_rows;
open my $list_fh, '<', $list_abs or die "Cannot read '$input_list': $!\n";
my $line_no = 0;
while (my $line = <$list_fh>) {
    ++$line_no;
    next if $line =~ /^\s*(?:#|$)/;
    chomp $line;
    my @fields = split /\s+/, $line;
    die "Invalid RCD row at line $line_no (expected: PDB start end chain sequence length CDR): $line\n"
        unless @fields >= 7;
    my ($pdb_field, $start, $end, $chain, $sequence, $length, $cdr) = @fields[0..6];
    die "Invalid start/end indices at line $line_no: $line\n"
        unless $start =~ /^-?\d+$/ && $end =~ /^-?\d+$/ && $end >= $start;
    die "Sequence length mismatch at line $line_no: declared $length, sequence has ".length($sequence)." residues.\n"
        unless $length =~ /^\d+$/ && $length == length($sequence);
    die "Invalid chain ID at line $line_no: '$chain'\n" unless length($chain) == 1;

    my $pdb_path = File::Spec->file_name_is_absolute($pdb_field)
        ? $pdb_field : File::Spec->catfile($list_dir, $pdb_field);
    $pdb_path = File::Spec->catfile(File::Spec->rel2abs('.'), $pdb_field)
        if !-f $pdb_path && -f $pdb_field;
    die "Structure from line $line_no not found: $pdb_field\n" unless -f $pdb_path;
    my $pdb_abs = abs_path($pdb_path);
    my $interval_length = $end - $start + 1;
    my $shift = $length - $interval_length;
    my $new_end = $start + $length - 1;
    my $out_path = $pdb_abs;
    $out_path =~ s/_renum\.pdb\z/_design.pdb/i
        or $out_path =~ s/\.pdb\z/_design.pdb/i
        or die "Input structure must have a .pdb extension: $pdb_field\n";

    my $key = join("\0", $pdb_abs, $start, $end, $chain, $sequence, $length, $cdr);
    $design_rows{$key} = {
        pdb => $pdb_abs, out => $out_path, pdb_field => File::Spec->abs2rel($out_path, $list_dir),
        start => $start, end => $end, new_end => $new_end, chain => $chain,
        sequence => $sequence, length => $length, cdr => $cdr, shift => $shift,
        line_no => $line_no,
    };
}
close $list_fh;
die "No RCD input rows found in '$input_list'.\n" unless %design_rows;

my %by_structure;
push @{$by_structure{$_->{pdb}}}, $_ for values %design_rows;
for my $pdb (sort keys %by_structure) {
    my @rows = @{$by_structure{$pdb}};
    my %seen_region;
    for my $row (@rows) {
        my $region = join(':', @{$row}{qw(chain start end)});
        die "Overlapping or duplicate design interval for '$pdb' at line $row->{line_no}.\n"
            if $seen_region{$region}++;
    }

    my $out = $rows[0]{out};
    open my $in_fh, '<', $pdb or die "Cannot read '$pdb': $!\n";
    open my $out_fh, '>', $out or die "Cannot write '$out': $!\n";
    while (my $line = <$in_fh>) {
        unless ($line =~ /^(?:ATOM  |HETATM)/) { print {$out_fh} $line; next; }
        my $chain = substr($line, 21, 1);
        my $resnum = substr($line, 22, 4);
        $resnum =~ s/\s//g;
        next unless $resnum =~ /^-?\d+$/;
        my $number = 0 + $resnum;
        my @matching = grep { $_->{chain} eq $chain && $number >= $_->{start} && $number <= $_->{end} } @rows;
        next if @matching;
        for my $row (@rows) {
            next unless $row->{chain} eq $chain && $number > $row->{end};
            my $new_number = $number + $row->{shift};
            die "Renumbered residue $new_number is outside the PDB residue-number field in '$pdb'.\n"
                if $new_number < -999 || $new_number > 9999;
            substr($line, 22, 4, sprintf('%4d', $new_number));
        }
        print {$out_fh} $line;
    }
    close $in_fh;
    close $out_fh or die "Error closing '$out': $!\n";
    print "Wrote ", basename($out), "\n";
}

my %cdr_rows;
for my $row (values %design_rows) {
    push @{$cdr_rows{$row->{cdr}}}, sprintf("%-20s %4d %4d %1s %20s %2d %s\n",
        $row->{pdb_field}, $row->{start}, $row->{new_end}, $row->{chain},
        $row->{sequence}, $row->{length}, $row->{cdr});
}
for my $cdr (sort keys %cdr_rows) {
    my $path = File::Spec->catfile($list_dir, "AbRaCD_design_${cdr}.txt");
    open my $out_fh, '>', $path or die "Cannot write '$path': $!\n";
    print {$out_fh} sort @{$cdr_rows{$cdr}};
    close $out_fh or die "Error closing '$path': $!\n";
    print "Wrote ", basename($path), "\n";
}

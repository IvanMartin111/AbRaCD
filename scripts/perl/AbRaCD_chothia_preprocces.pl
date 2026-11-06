#!/usr/bin/perl

# Renumber Chothia antibody structures and prepare one RCD input per CDR.
use strict;
use warnings;
use File::Basename qw(dirname basename);
use File::Spec;
use Cwd qw(abs_path);

my %aa = (
    ALA=>'A', CYS=>'C', CYX=>'C', ASP=>'D', GLU=>'E', PHE=>'F', GLY=>'G',
    HIS=>'H', HSD=>'H', HSE=>'H', HSP=>'H', ILE=>'I', LYS=>'K', LEU=>'L',
    MET=>'M', ASN=>'N', PRO=>'P', GLN=>'Q', ARG=>'R', SER=>'S', THR=>'T',
    VAL=>'V', TRP=>'W', TYR=>'Y',
);
my %cdr = (
    L1 => ['L',24,34], L2 => ['L',50,56], L3 => ['L',89,97],
    H1 => ['H',26,32], H2 => ['H',52,56], H3 => ['H',93,102],
);

sub usage {
    print <<'USAGE';
AbRaCD Chothia preprocessing
============================

USAGE
    print "  AbRaCD_chothia_preprocces.pl [OPTIONS] [PDB|PDB_LIST.txt]\n\n";

    print "INPUT\n";
    print "  No input       Process every original *.pdb in the current directory.\n";
    print "  PDB            Process only that structure.\n";
    print "  PDB_LIST.txt   Process only the PDB files listed in the text file.\n";
    print "                 Missing files are reported and skipped.\n\n";

    print "CHAIN SELECTION WITHOUT --purge\n";
    print "  Accepted: H and L; H/L plus one antigen chain; or one of H/L plus one antigen.\n";
    print "  Antigen and other unselected chains are removed. H and L are identified by\n";
    print "  their chain IDs.\n";

    print "OPTIONS\n";
    print "  --purge FILE   Use FILE to select the PDBs and read their source H/L chain IDs.\n";
    print "                 Keep only those chains and rename them H and L in processed PDBs.\n\n";

    print "PURGE FILE FORMAT\n";
    print "  One whitespace-separated row per structure:  PDB-path  H-chain  L-chain\n";
    print "  With --purge, only PDBs listed in this file are processed. Relative PDB paths\n";
    print "  are resolved from the directory containing the purge file.\n";
    print "  Use '-' for an absent chain (for example, a nanobody). Lines starting with #\n";
    print "  and blank lines are ignored. PDB-name may include or omit the .pdb suffix.\n\n";
    print "  Example:\n";
    print "    pdb1 A B\n";
    print "    pdb2 B Z\n";
    print "    nanobody.pdb A -\n\n";
    print "  Example command:\n";
    print "    AbRaCD_chothia_preprocces.pl --purge purge.txt\n\n";

    print "PROCESSING AND OUTPUTS\n";
    print "  All six CDR inputs are prepared. Only selected H/L chains are retained.\n";
    print "  Residues after original Chothia H113 or L109 are removed before renumbering.\n";
    print "  Renumbered PDBs and AbRaCD_input_{CDR}.txt files are written beside the\n";
    print "  input structures. PyRosetta is required for open_and_close_batch_par.py.\n";
}

if (@ARGV && ($ARGV[0] eq '-h' || $ARGV[0] eq '--help')) { usage(); exit 0; }
my $purge_file;
if (@ARGV && $ARGV[0] eq '--purge') {
    shift @ARGV;
    die "--purge requires a purge.txt file.\n" unless @ARGV;
    $purge_file = shift @ARGV;
}
die "Usage: AbRaCD_chothia_preprocces.pl [--purge purge.txt] [structure.pdb|structures.txt]\n" if @ARGV > 1;
my $input = @ARGV ? shift @ARGV : undef;
$input = $purge_file if defined $purge_file;

my @pdbs;
if (defined $purge_file) {
    open my $purge_list, '<', $purge_file or die "Cannot read purge file '".basename($purge_file)."': $!\n";
    my $base = dirname(abs_path($purge_file));
    my @missing;
    while (<$purge_list>) {
        chomp; s/^\s+|\s+$//g; next if !length || /^#/;
        my ($name) = split /\s+/;
        my $path = File::Spec->file_name_is_absolute($name) ? $name : File::Spec->catfile($base, $name);
        if (-f $path) { push @pdbs, abs_path($path); }
        else { push @missing, basename($path); }
    }
    close $purge_list;
    if (@missing) {
        print STDERR "WARNING: the following PDB(s) from '$purge_file' were not found and will be skipped:\n";
        print STDERR "  - $_\n" for @missing;
    }
} elsif (!defined $input) {
    @pdbs = map { abs_path($_) } grep { -f $_ && $_ !~ /_renum\.pdb$/i } glob('*.pdb');
    die "No original .pdb files found in the current directory.\n" unless @pdbs;
} elsif ($input =~ /\.txt$/i) {
    open my $list, '<', $input or die "Cannot read list '".basename($input)."': $!\n";
    my $base = dirname(abs_path($input));
    my @missing;
    while (<$list>) {
        chomp; s/^\s+|\s+$//g; next if !length || /^#/;
        my $path = File::Spec->file_name_is_absolute($_) ? $_ : File::Spec->catfile($base, $_);
        if (-f $path) { push @pdbs, abs_path($path); }
        else { push @missing, basename($path); }
    }
    close $list;
    if (@missing) {
        print STDERR "WARNING: the following PDB(s) from '$input' were not found and will be skipped:\n";
        print STDERR "  - $_\n" for @missing;
    }
} else { push @pdbs, $input; }
die "No PDB structures found in '$input'.\n" if defined($input) && !@pdbs;

# Resolve the chain mapping for every structure before making any output.
my %purge_rows;
if (defined $purge_file) {
    open my $purge, '<', $purge_file or die "Cannot read purge file '".basename($purge_file)."': $!\n";
    while (<$purge>) {
        chomp; s/^\s+|\s+$//g; next if !length || /^#/;
        my ($name,$heavy,$light) = split /\s+/;
        die "Invalid purge.txt row (expected: PDB H-chain L-chain): $_\n" unless defined $light;
        my $base=basename($name); my $stem=$base; $stem =~ s/\.pdb$//i;
        $purge_rows{$base}=[$heavy,$light]; $purge_rows{$stem}=[$heavy,$light];
    }
    close $purge;
}
my %chain_map;
for my $pdb (@pdbs) {
    my $pdb_name=basename($pdb);
    open my $fh, '<', $pdb or die "Cannot read PDB '$pdb_name': $!\n";
    my %chains;
    while (my $line=<$fh>) { next unless $line =~ /^ATOM\s/; my $c=substr($line,21,1); $chains{$c}=1; }
    close $fh;
    my @chains=sort keys %chains;
    my %mapping;
    if (defined $purge_file) {
        my $row=$purge_rows{basename($pdb)};
        my $stem=basename($pdb); $stem =~ s/\.pdb$//i;
        $row ||= $purge_rows{$stem};
        die "No --purge mapping for '$pdb_name' in '".basename($purge_file)."'.\n" unless $row;
        my ($heavy,$light)=@$row;
        $heavy='-' if !defined($heavy) || $heavy =~ /^(?:NA|NONE)$/i;
        $light='-' if !defined($light) || $light =~ /^(?:NA|NONE)$/i;
        die "Purge mapping for '$pdb_name' must name distinct H/L source chains (or use - for a missing chain).\n"
            if $heavy ne '-' && $light ne '-' && $heavy eq $light;
        die "Purge mapping for '$pdb_name' names a chain that is absent from the PDB.\n"
            if ($heavy ne '-' && !$chains{$heavy}) || ($light ne '-' && !$chains{$light});
        $mapping{$heavy}='H' if $heavy ne '-';
        $mapping{$light}='L' if $light ne '-';
        die "Purge mapping for '$pdb_name' must retain at least one chain.\n" unless %mapping;
        print "WARNING: '$pdb_name' has no heavy chain selected by --purge.\n" if $heavy eq '-';
        print "WARNING: '$pdb_name' has no light chain selected by --purge.\n" if $light eq '-';
    } else {
        my $has_h=$chains{'H'}; my $has_l=$chains{'L'};
        my $valid = (@chains==2 && $has_h && $has_l)
                 || (@chains==2 && ($has_h || $has_l))
                 || (@chains==3 && $has_h && $has_l);
        die "ERROR: '$pdb_name' has chains [".join(' ',@chains)."]. Expected H and L, H/L plus one antigen, or H plus one antigen. Use --purge purge.txt (PDB H-chain L-chain) to select chains. No outputs were created.\n"
            unless $valid;
        $mapping{H}='H' if $has_h;
        $mapping{L}='L' if $has_l;
        print "WARNING: '$pdb_name' has no H chain; only the L chain will be processed.\n" if !$has_h && $has_l;
        print "WARNING: '$pdb_name' has no L chain; only the H chain will be processed.\n" if $has_h && !$has_l;
    }
    $chain_map{$pdb}=\%mapping;
    my @removed=grep { !exists $mapping{$_} } @chains;
    if (@removed) {
        print "WARNING: '$pdb_name': removed ",scalar(@removed)," chain(s): ",join(', ',@removed),".\n";
    }
}

my %inputs = map { $_ => [] } keys %cdr;
my %gap_inputs = map { $_ => [] } keys %cdr;
my @oc_names;
my %oc_dirs;
my $failed_gap = 0;

for my $pdb (@pdbs) {
    my $pdb_name=basename($pdb);
    open my $fh, '<', $pdb or die "Cannot read PDB '$pdb_name': $!\n";
    my (%is_h, %is_l, %residues);
    $is_h{H}=1 if grep { $chain_map{$pdb}{$_} eq 'H' } keys %{$chain_map{$pdb}};
    $is_l{L}=1 if grep { $chain_map{$pdb}{$_} eq 'L' } keys %{$chain_map{$pdb}};
    while (my $line = <$fh>) {
        next unless $line =~ /^ATOM\s/;
        my $source_chain = substr($line,21,1);
        next unless exists $chain_map{$pdb}{$source_chain};
        my $chain = $chain_map{$pdb}{$source_chain};
        my $num = substr($line,22,4); $num =~ s/\s//g;
        my ($original_number) = $num =~ /^(-?\d+)/;
        my $max_number = $chain eq 'H' ? 113 : 109;
        next if defined($original_number) && $original_number > $max_number;
        my $icode = substr($line,26,1); my $resname = substr($line,17,3); $resname =~ s/\s//g;
        my $key = join(':',$chain,$num,$icode);
        if (!exists $residues{$key}) {
            $residues{$key} = { chain=>$chain, num=>$num, icode=>$icode, res=>$resname, lines=>[], ca=>undef };
        }
        my $atom = substr($line,12,4); my $element = substr($line,76,2); $element =~ s/\s//g;
        next if $element eq 'H' || $atom =~ /^\s*H/;
        push @{$residues{$key}{lines}}, $line;
        if ($atom =~ /^\s*CA\s*$/) {
            $residues{$key}{ca} = [map { 0+substr($line,$_,8) } (30,38,46)];
        }
    }
    close $fh;
    my @ordered = sort {
        ($a->{chain} cmp $b->{chain}) ||
        ((($a->{num}=~/^-?\d+$/ && $b->{num}=~/^-?\d+$/) ? ($a->{num}<=>$b->{num}) : ($a->{num} cmp $b->{num}))) ||
        ($a->{icode} cmp $b->{icode})
    } values %residues;

    my %newnum; my %loop_members; my @gaps; my %chain_count;
    my %last;
    for my $r (@ordered) {
        my $chain=$r->{chain};
        $newnum{join(':',$chain,$r->{num},$r->{icode})} = ++$chain_count{$chain};
        if (exists $last{$chain}) {
            my $prev=$last{$chain};
            my ($pn) = $prev->{num} =~ /^(-?\d+)/; my ($cn) = $r->{num} =~ /^(-?\d+)/;
            my $missing = defined($pn) && defined($cn) && $cn-$pn>1 ? $cn-$pn-1 : 0;
            my $dist;
            if ($prev->{ca} && $r->{ca}) {
                $dist=sqrt(($prev->{ca}[0]-$r->{ca}[0])**2+($prev->{ca}[1]-$r->{ca}[1])**2+($prev->{ca}[2]-$r->{ca}[2])**2);
            }
            # A Chothia numbering jump alone is not a structural gap. Require
            # both missing residue numbers and a long CA-to-CA distance.
            if ($missing && defined($dist) && $dist>4.2) {
                my @at_cdr;
                for my $name (keys %cdr) {
                    my ($c,$lo,$hi)=@{$cdr{$name}};
                    next unless (($c eq 'H' && $is_h{$chain}) || ($c eq 'L' && $is_l{$chain}));
                    my ($pnum)=$prev->{num}=~/^(\d+)/; my ($cnum)=$r->{num}=~/^(\d+)/;
                    next unless defined($pnum) && defined($cnum);
                    push @at_cdr,$name if ($pnum >= $lo && $pnum <= $hi) || ($cnum >= $lo && $cnum <= $hi);
                }
                my $new_prev=$newnum{join(':',$chain,$prev->{num},$prev->{icode})};
                my $new_curr=$newnum{join(':',$chain,$r->{num},$r->{icode})};
                push @gaps, [$chain,$prev->{num}.$prev->{icode},$r->{num}.$r->{icode},$missing,$new_prev,$new_curr,\@at_cdr];
            }
        }
        $last{$chain}=$r;
        for my $name (keys %cdr) {
            my ($c,$lo,$hi)=@{$cdr{$name}};
            next unless (($c eq 'H' && $is_h{$chain}) || ($c eq 'L' && $is_l{$chain}));
            my ($n)=$r->{num}=~/^(\d+)/;
            push @{$loop_members{$name}},$r if defined($n) && $n >= $lo && $n <= $hi;
        }
    }

    my %critical_cdr;
    for my $g (@gaps) {
        my ($chain,$from,$to,$missing,$new_from,$new_to,$inside)=@$g;
        my @predict=@$inside;
        if (@predict) {
            print "ERROR: $pdb_name chain $chain has a CRITICAL STRUCTURAL GAP in predicted CDR ",join(',',@predict)," ($from -$to ). PLEASE RECONSIDER; use the antibody-engineering preparation script.\n";
            $critical_cdr{$_}=1 for @predict;
            $failed_gap=1;
        } else {
            printf "WARNING: %s chain %s structural gap of %d residue(s): original %s-%s; after renumbering %s-%s. Continuing.\n", $pdb_name,$chain,$missing,$from,$to,$new_from,$new_to;
        }
    }

    my $out = $pdb; $out =~ s/\.pdb$/_renum.pdb/i;
    die "Input '$pdb_name' must have a .pdb extension.\n" if $out eq $pdb;
    open my $outfh, '>', $out or die "Cannot write '".basename($out)."': $!\n";
    print {$outfh} "REMARK   Renumbered consecutively by AbRaCD_chothia_preprocces.pl\n";
    my %atom_serial;
    for my $r (@ordered) {
        my $key=join(':',$r->{chain},$r->{num},$r->{icode});
        for my $line (@{$r->{lines}}) {
            substr($line,6,5,sprintf('%5d',++$atom_serial{$r->{chain}}));
            substr($line,21,1,$r->{chain});
            substr($line,22,4,sprintf('%4d',$newnum{$key})); substr($line,26,1,' ');
            print {$outfh} $line;
        }
    }
    print {$outfh} "END\n"; close $outfh;
    my $absout=abs_path($out);
    push @oc_names,$absout;
    $oc_dirs{dirname($absout)}=1;

    for my $name (keys %cdr) {
        next unless @{$loop_members{$name} || []};
        my ($chain)=@{$cdr{$name}};
        my ($chain_id)=map { $_->{chain} } @{$loop_members{$name}};
        my $first=$loop_members{$name}[0]; my $last=$loop_members{$name}[-1];
        my ($start,$end)=map { $newnum{join(':',$_->{chain},$_->{num},$_->{icode})} } ($first,$last);
        my $seq=join('',map { $aa{$_->{res}} // 'X' } @{$loop_members{$name}});
        if ($name eq 'L3' && (length($seq) == 9 || length($seq) == 10)) {
            my $cis_proline_pos = length($seq) - 3;
            substr($seq,$cis_proline_pos,1,'p') if substr($seq,$cis_proline_pos,1) eq 'P';
        }
        my $chain_field=$chain_id eq ' ' ? ' ' : $chain_id;
        # Keep the whitespace-separated RCD fields in fixed-width columns for
        # readability. The 20-character sequence field is right-aligned so
        # shorter loops line up at their C-terminal end.
        my $input_row = sprintf("%-20s %4d %4d %1s %20s %2d %s",basename($absout),$start,$end,$chain_field,$seq,length($seq),$name);
        if ($critical_cdr{$name}) {
            push @{$gap_inputs{$name}}, $input_row;
        } else {
            push @{$inputs{$name}}, $input_row;
        }
    }
}

# RCD lists are written beside the supplied input list, or beside the single PDB.
my $list_dir = defined($input) && $input =~ /\.txt$/i ? dirname(abs_path($input)) : dirname(abs_path($pdbs[0]));
for my $name (sort keys %cdr) {
    my $path=File::Spec->catfile($list_dir,"AbRaCD_input_${name}.txt");
    open my $fh,'>',$path or die "Cannot write '".basename($path)."': $!\n";
    print {$fh} "$_\n" for @{$inputs{$name}};
    close $fh;

    next unless @{$gap_inputs{$name}};
    my $gap_path=File::Spec->catfile($list_dir,"AbRaCD_GAP_${name}.txt");
    open my $gap_fh,'>',$gap_path or die "Cannot write '".basename($gap_path)."': $!\n";
    print {$gap_fh} "$_\n" for @{$gap_inputs{$name}};
    close $gap_fh;
}

# PyRosetta batch round-trip. Keep only heavy ATOM and TER records afterward.
if (!@oc_names) {
    print "No structures were renumbered; no open_and_close batch was run.\n";
    exit($failed_gap ? 2 : 1);
}
my @dirs=sort keys %oc_dirs;
die "All structures must be in one directory for the batch open_and_close output.\n" if @dirs>1;
my $structure_list=File::Spec->catfile($dirs[0],'structures.txt');
open my $listfh,'>',$structure_list or die "Cannot write 'structures.txt': $!\n";
print {$listfh} basename($_),"\n" for @oc_names; close $listfh;
my $old_dir=File::Spec->rel2abs('.');
chdir $dirs[0] or die "Cannot enter the structure directory for open_and_close: $!\n";
my $status=system('open_and_close_batch_par.py', basename($structure_list), '_oc');
chdir $old_dir or die "Cannot return to the starting directory: $!\n";
if ($status != 0) { die "open_and_close_batch_par.py failed; temporary _oc files and structures.txt were kept for diagnosis.\n"; }
for my $path (@oc_names) {
    my $oc=$path; $oc =~ s/\.pdb$/_oc.pdb/i;
    open my $in,'<',$oc or die "Expected open_and_close output '".basename($oc)."' was not created.\n";
    my $tmp="$path.tmp"; open my $clean,'>',$tmp or die "Cannot write '".basename($tmp)."': $!\n";
    while (my $line=<$in>) {
        next unless $line =~ /^(?:ATOM\s|TER\s|TER\n|TER$)/;
        next if $line =~ /^ATOM\s/ && (substr($line,76,2) =~ /^\s*H/i || substr($line,12,4) =~ /^\s*H/);
        print {$clean} $line;
    }
    close $in; close $clean;
    rename $tmp,$path or die "Cannot replace '".basename($path)."': $!\n";
    unlink $oc or warn "Could not remove temporary '".basename($oc)."': $!\n";
}
unlink $structure_list or warn "Could not remove 'structures.txt': $!\n";
print "Prepared ",scalar(@pdbs)," structure(s). Renumbered PDBs and per-CDR RCD lists are beside the source structure(s).\n";
exit($failed_gap ? 2 : 0);

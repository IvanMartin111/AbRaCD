# AbRaCD command-line manual

This manual describes the command-line interface implemented by `scripts/perl/abracd.pl` and the input files it consumes. Run commands from the working directory where you want the generated files, result folders, and log. The launcher can be checked with `abracd.pl --help`.

## 1. Requirements and setup checks

Use the installation steps in [Readme.md](Readme.md). Before a run, confirm that:

- `abracd.pl` is on `PATH` and Perl can execute it.
- Python can import `pyrosetta`. 
- MPI runs also require `mpi4py` and `mpirun`.
- `data/Maps/CDRS.bin`, `data/Maps/korp6Dv1.bin`, and `data/Maps/dunbrack.bin` exist.
- An RCD executable is on `PATH`. Serial selection looks for `rcd_gcc` or `rcd`; MPI runs look for `rcd_mpi_gcc`, `rcd_mpi_intel`, or `rcd_mpi` as well as `mpirun`. If only an MPI executable is available, the launcher can fall back to it.
- Rosetta refinement helpers are available through the project installation. `--resROS 0` skips Rosetta refinement.

The launcher checks several dependencies and input files before starting. Its log records the invocation and executed commands.

## 2. Main command

```text
abracd.pl [INPUT] [--raw | --model] [--purge] --CDRs CDR [CDR ... | ALL] [OPTIONS]
```

`INPUT` is optional except where a mode requires it. `--CDRs` is required. Options are case-sensitive. The mode flag must immediately follow the optional input (or appear first when there is no input). The selected CDRs run in the order given. `ALL` expands to `H1 H2 H3 L1 L2 L3` and must be used alone.

### Input modes

| Mode | Input and behavior |
|---|---|
| Default preprocessing | `INPUT` may be one `.pdb`, or a `.txt` list of PDB paths. With no input, all original `*.pdb` files in the current directory are processed. Chothia preprocessing runs before the requested CDRs. |
| `--purge` | Provide the purge table as `INPUT`; it selects structures and source H/L chain IDs. Only valid in default preprocessing mode. |
| `--raw` | Skips Chothia preprocessing. For one CDR, provide that CDR's RCD-format input list as `INPUT`. For multiple CDRs or `ALL`, omit `INPUT`; the launcher reads `AbRaCD_input_<CDR>.txt` files from the current directory. |
| `--model` | Design/modeling mode without native-structure RMSD. Provide one RCD-format `.txt` input and exactly one CDR. The launcher prepares design PDBs, runs RCD and refinement, and removes native-dependent score columns. Cannot be combined with `--raw` or `--purge`. |

### Options

| Option | Value | Meaning and default |
|---|---|---|
| `--CDRs` | one or more of `H1 H2 H3 L1 L2 L3`, or `ALL` | Required. CDRs are processed in the order listed. No duplicates; `ALL` must be the sole value. |
| `--purge` | none | Treat the positional input as a purge table. Preprocessing retains the specified chains and renames them to H and L. |
| `--raw` | none | Use prepared RCD input files and skip Chothia preprocessing. See the mode table for the positional-input rule. |
| `--model` | none | Prepare and run design/modeling workflow. Requires one RCD list and one CDR. |
| `--OUT` | `NAME` | Output-name prefix, optional. Allowed characters: letters, digits, `.`, `_`, `-`. The log becomes `NAME.log`; generated RCD output names include the prefix. Default: no prefix. |
| `--NP` | positive integer | Request MPI RCD with this many processes. Without this option, serial RCD is preferred with one process; if only an MPI RCD executable is available, the launcher falls back to MPI with the default process count of `1`. MPI and `mpirun` must be installed for MPI execution. |
| `--NLOOPS` | positive integer | Number of loop samples requested from RCD. Defaults: H3 `50000`; other CDRs `20000`. |
| `--LOCOBEST` | positive integer | RCD `loco_best` value. Defaults: H3 `1000`; other CDRs `2000`. |
| `--repRCD` | positive integer | Number of RCD sampling runs per CDR. Default: `1`. |
| `--resROS` | positive integer or `0` | Number of Rosetta refinement runs per RCD run. Default: `1`; `0` disables refinement. |

Each option may be specified at most once. Integer parameters must be positive, except `--resROS` accepts zero. `--purge` does not accept a value in the launcher: it uses the same file supplied as `INPUT`.

### Examples

```bash
# One Chothia-numbered antibody structure, one CDR
abracd.pl antibody.pdb --CDRs H3

# A text list of PDBs, several selected CDRs
abracd.pl antibodies.txt --CDRs H1 H3 L1

# Process all six CDRs; default processing of PDBs in the current directory
abracd.pl --CDRs ALL

# Explicit source chain IDs for each structure
abracd.pl purge.txt --purge --CDRs H3 --NP 8 --OUT purified_h3

# Use a pre-existing RCD list for one CDR
abracd.pl AbRaCD_input_H3.txt --raw --CDRs H3 --NLOOPS 30000

# Use pre-existing lists for multiple CDRs
abracd.pl --raw --CDRs H1 H3

# Design/modeling mode for one loop
abracd.pl AbRaCD_GAP_H3.txt --model --CDRs H3 --NP 8 --OUT design_h3

# Repeated sampling and no Rosetta refinement
abracd.pl antibody.pdb --CDRs H3 --repRCD 3 --resROS 0
```

Invalid combinations include `--raw` and `--model` together, `--purge` with either mode, a positional input for `--raw` with multiple CDRs, and multiple CDRs with `--model`.

## 3. File formats

### PDB list

A `.txt` structure list contains one PDB path per line. Relative paths are resolved from the directory containing the list. Blank lines are ignored. In default preprocessing mode, use either one PDB directly, or this list.

```text
structures/antibody_a.pdb
structures/antibody_b.pdb
```

When no input is provided, the preprocessor scans the current directory for original `.pdb` files (excluding files ending in `_renum.pdb`).

### Purge table

A purge table contains one whitespace-separated row per structure:

```text
PDB_path  H_chain_ID  L_chain_ID
```

Relative PDB paths are resolved from the purge file's directory. Use `-` for a missing chain, such as a nanobody. Blank lines and lines beginning with `#` are ignored. The structure name may be written with or without `.pdb` in the mapping lookup. The table defines which source chains are retained and mapped to H/L.

```text
# structure heavy-chain light-chain
complex_a.pdb A B
complex_b B Z
nanobody.pdb H -
```

Run it through the main workflow as `abracd.pl purge.txt --purge --CDRs H3`, or invoke preprocessing directly as described below.

### RCD input list

Each non-comment row has at least seven whitespace-separated fields:

```text
PDB_file  start  end  chain  sequence  length  CDR
```

`start` and `end` are inclusive residue indices; `end` must be greater than or equal to `start`. The sequence length must equal `length`. Chain is a one-character chain ID. Relative PDB paths are resolved from the RCD list directory (the design-preparation helper also accepts a matching path in the current directory).

```text
antibody_renum.pdb  97  109  H  ARWGGDGFYAMDY  13  H3
```

Comments and blank lines are ignored by the design preparation helper. Standard generated files are named `AbRaCD_input_<CDR>.txt`. For `--model`, the input is expected to describe the loop/design interval in this same RCD row format; the helper writes `AbRaCD_design_<CDR>.txt` and `<source>_design.pdb` beside the input structures/list.

## 4. Direct helper commands

These scripts can be called independently when only a preparation or utility step is needed.

### Chothia preprocessing

```text
AbRaCD_chothia_preprocces.pl [--purge PURGE.txt] [PDB | PDB_LIST.txt]
```

Examples:

```bash
AbRaCD_chothia_preprocces.pl antibody.pdb
AbRaCD_chothia_preprocces.pl structures.txt
AbRaCD_chothia_preprocces.pl --purge purge.txt
AbRaCD_chothia_preprocces.pl --help
```

With no input, it processes original PDBs in the current directory. It generates renumbered structures and RCD input files for all six CDRs beside the input/list. Without a purge file, accepted chain layouts are H and L, H/L plus one antigen, or one H/L chain plus one antigen. Other unselected chains are removed. Residues after original Chothia H113 or L109 are removed before renumbering. PyRosetta is used by the structure-loading helper.

### Engineering preparation

```text
ab_engineering_prep.pl RCD_INPUT.txt
```

Example:

```bash
ab_engineering_prep.pl AbRaCD_GAP_H3.txt
```

For each input row it removes the specified loop residues from that chain and shifts subsequent residue numbering to fit the requested sequence length. It writes `<source>_design.pdb` and `AbRaCD_design_<CDR>.txt`. The main `--model` workflow calls this helper automatically.

### PDB load/save utility

```text
open_and_close_batch_par.py INPUT OUTPUT [LIGAND_PARAMS] [--par N]
```

`INPUT` is one PDB or a `.txt` list of PDB paths; `OUTPUT` is the output path for one PDB, or a suffix for each listed PDB (for example `_oc`). Optional `LIGAND_PARAMS` is a comma-separated list of ligand `.params` files. `--par N` uses N MPI processes; N defaults to 1 (sequential) and must be at least 1.

```bash
open_and_close_batch_par.py antibody.pdb antibody_oc.pdb
open_and_close_batch_par.py structures.txt _oc --par 8
open_and_close_batch_par.py complex.pdb complex_oc.pdb ligand.params
```

The utility loads and writes structures through PyRosetta; MPI mode additionally needs `mpi4py` and `mpiexec`.

## 5. Outputs and logs

- The launcher writes a detailed run log in the current directory: `AbRaCD.log`, or `<OUT>.log` when `--OUT NAME` is provided.
- Preprocessing writes renumbered PDBs and per-CDR input lists beside the source/list directory.
- RCD creates a result folder for each CDR/run, with its name derived from `--OUT` and the CDR/repetition suffixes.
- Rosetta refinement writes refined structures and score tables in the working layout used by `refine_loops.pl`.
- Modeling mode removes result sections before `RUN SPECIFICATIONS` and strips native-dependent RMSD columns from Rosetta tables, because no native reference is available.

Check the log for the exact executed commands, stage names, and failure details. Keep the launch working directory consistent when you need predictable output placement.

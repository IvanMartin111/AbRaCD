#!/home/ivan/miniconda3/envs/igfold/bin/python

##!/usr/bin/python3

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Load PDB files with PyRosetta, sequentially or in parallel with MPI.")
    parser.add_argument("input", help="PDB file or .txt list of PDB files")
    parser.add_argument("output", help="output PDB path, or suffix when input is a .txt list")
    parser.add_argument("params", nargs="?", default="", help="comma-separated ligand .params files")
    parser.add_argument("--par", type=int, default=1, metavar="N", help="number of MPI processes (default: 1, sequential)")
    parser.add_argument("--mpi-worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.par < 1:
        parser.error("--par must be at least 1")
    return args


def get_jobs(input_path, output_arg):
    if input_path.suffix.lower() == ".txt":
        names = [line.strip() for line in input_path.read_text().splitlines() if line.strip()]
        if not names:
            raise ValueError("Input list is empty: %s" % input_path)
        jobs = []
        for name in names:
            pdb_path = Path(name)
            if not pdb_path.is_absolute():
                pdb_path = input_path.parent / pdb_path
            output_path = pdb_path.with_name(pdb_path.stem + output_arg + pdb_path.suffix)
            jobs.append((str(pdb_path), str(output_path)))
        return jobs
    return [(str(input_path), output_arg)]


def process_pdb(pdb_path, output_path, params):
    from pyrosetta import pose_from_file
    if params:
        from pyrosetta import Vector1, generate_nonstandard_residue_set
        restype = generate_nonstandard_residue_set(Vector1(params.split(",")))
        pose = pose_from_file(restype, pdb_path)
    else:
        pose = pose_from_file(pdb_path)
    pose.dump_pdb(output_path)
    return pose.total_residue()


def initialize_pyrosetta():
    from pyrosetta import init
    init(extra_options="-ex2 -rebuild_disulf false -detect_disulf false -ignore_zero_occupancy false -allow_omega_move false -ignore_unrecognized_res -out:levels core:error -out:levels protocols:error")


def worker(args):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    jobs = get_jobs(Path(args.input), args.output)
    initialize_pyrosetta()
    local_results = []
    batch_start = time.perf_counter()
    for index in range(rank, len(jobs), size):
        pdb_path, output_path = jobs[index]
        start = time.perf_counter()
        print("[rank %d] Loading %s" % (rank, pdb_path), flush=True)
        try:
            residues = process_pdb(pdb_path, output_path, args.params)
            elapsed = time.perf_counter() - start
            print("[rank %d] %s -> %s (%d residues): %.3f s" % (rank, pdb_path, output_path, residues, elapsed), flush=True)
            local_results.append((index, pdb_path, elapsed, None))
        except Exception as exc:
            elapsed = time.perf_counter() - start
            print("[rank %d] ERROR %s after %.3f s: %s" % (rank, pdb_path, elapsed, exc), file=sys.stderr, flush=True)
            local_results.append((index, pdb_path, elapsed, str(exc)))
    gathered_results = comm.gather(local_results, root=0)
    if rank == 0:
        all_results = [item for group in gathered_results for item in group]
        all_results.sort()
        total = sum(item[2] for item in all_results)
        average = total / len(all_results) if all_results else 0.0
        wall = time.perf_counter() - batch_start
        print("\nResumen: %d casos; suma de tiempos por caso %.3f s; promedio %.3f s/caso; tiempo de pared aprox. %.3f s" % (len(all_results), total, average, wall), flush=True)
        failures = sum(item[3] is not None for item in all_results)
        if failures:
            print("Casos con error: %d" % failures, file=sys.stderr, flush=True)
            return 1
    return 0


def main():
    args = parse_args()
    if args.mpi_worker:
        return worker(args)
    input_path = Path(args.input)
    if args.par == 1:
        jobs = get_jobs(input_path, args.output)
        initialize_pyrosetta()
        durations = []
        failures = 0
        total_start = time.perf_counter()
        for pdb_path, output_path in jobs:
            start = time.perf_counter()
            try:
                print("Loading %s" % pdb_path, flush=True)
                residues = process_pdb(pdb_path, output_path, args.params)
                elapsed = time.perf_counter() - start
                print("%s -> %s (%d residues): %.3f s" % (pdb_path, output_path, residues, elapsed), flush=True)
            except Exception as exc:
                elapsed = time.perf_counter() - start
                failures += 1
                print("ERROR %s after %.3f s: %s" % (pdb_path, elapsed, exc), file=sys.stderr, flush=True)
            durations.append(elapsed)
        wall = time.perf_counter() - total_start
        total = sum(durations)
        average = total / len(durations) if durations else 0.0
        print("\nResumen: %d casos; suma de tiempos por caso %.3f s; promedio %.3f s/caso; tiempo total %.3f s" % (len(durations), total, average, wall), flush=True)
        return 1 if failures else 0

    command = ["mpiexec", "-n", str(args.par), sys.executable, str(Path(__file__).resolve()), args.input, args.output]
    if args.params:
        command.append(args.params)
    command.extend(["--par", str(args.par), "--mpi-worker"])
    env = os.environ.copy()
    env.setdefault("OMPI_ALLOW_RUN_AS_ROOT", "1")
    env.setdefault("OMPI_ALLOW_RUN_AS_ROOT_CONFIRM", "1")
    return subprocess.call(command, env=env)


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3

import os
import argparse
import pandas as pd


COMPUTE_COVERAGE_TEMPLATE = """bowtie2 \
-x {database} \
{reads} \
--end-to-end \
--very-sensitive \
--phred33 \
--threads {threads} \
--seed 0 \
--time \
-k 2 \
--no-unal \
--no-discordant \
-X {fragment} \
2>> {log} | \
sambamba view -q --nthreads {threads} \
--compression-level 6 \
--format bam \
--compression-level 6 \
--sam-input /dev/stdin \
--output-filename /dev/stdout | \
sambamba sort -q --nthreads {threads} \
--memory-limit {memory_limit} \
--compression-level 0 \
--tmpdir {temp_dir} \
--out /dev/stdout /dev/stdin | \
jgi_summarize_bam_contig_depths \
--outputDepth {coverage} - \
2>> {log}"""


def parse_samples(samples_tsv):
    return pd.read_csv(samples_tsv, sep="\s+").set_index("id", drop=False)


def get_fqpath(sample_df, sample_id, col):
    return sample_df.loc[sample_id, [col]].dropna()[0]


class coverager:
    def __init__(
        self, reads, database, fragment, threads, memory_limit, temp_dir, coverage, log,
    ):
        if len(reads) == 2:
            self.reads = "-1 %s -2 %s" % (reads[0], reads[1])
        else:
            self.reads = "-U %s" % reads[0]

        self.threads = threads
        self.database = database
        self.memory_limit = memory_limit
        self.temp_dir = temp_dir
        self.coverage = coverage
        self.fragment = fragment
        self.log = log


def main():
    parser = argparse.ArgumentParser(description="metapi-profiling.py")
    parser.add_argument("-s", "--samples", type=str, required=True, help="samples.tsv")
    parser.add_argument(
        "-d", "--database", type=str, default=None, help="host index base path"
    )
    parser.add_argument(
        "-l", "--fragment", type=int, default=1200, help="align fragment length"
    )
    parser.add_argument(
        "-m", "--memory_limit", type=str, default="8G", help="memory limit"
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=8, help="bowtie2, samamba threads",
    )
    parser.add_argument("-o", "--out_dir", type=str, help="output directory")

    args = parser.parse_args()

    samples_df = parse_samples(args.samples)

    os.makedirs(args.out_dir, exist_ok=True)

    with open("coverage.sh", "w") as oh:
        for sample_id in samples_df.index:
            reads = []
            reads.append(get_fqpath(samples_df, sample_id, "fq1"))
            if "fq2" in samples_df.columns:
                reads.append(get_fqpath(samples_df, sample_id, "fq2"))

            outdir = os.path.join(args.out_dir, sample_id)
            os.makedirs(outdir, exist_ok=True)
            temp_dir = os.path.join(outdir, "temp")
            coverage = os.path.join(outdir, sample_id + ".jgi.coverage.txt")
            log = os.path.join(outdir, sample_id + ".jgi.coverage.log")

            coverage_cmd = COMPUTE_COVERAGE_TEMPLATE.format_map(
                vars(
                    coverager(
                        reads,
                        args.database,
                        args.fragment,
                        args.threads,
                        args.memory_limit,
                        temp_dir,
                        coverage,
                        log,
                    )
                )
            )
            oh.write(coverage_cmd + "\n")

    os.chmod("coverage.sh", 0o755)


if __name__ == "__main__":
    main()

#!/usr/bin/env bash

set -euo pipefail

die() {
    echo $1 1>&2
    exit 1
}

nargs=3
if [ $# -ne ${nargs} ]; then
    die "Usage: $0 input_file head_lines tail_lines"
fi
input_file=$1
head_lines=$2
tail_lines=$3

output_file="chunk.txt"

workdir=$(mktemp -d)
head_output="${workdir}/temp.txt"
head --lines ${head_lines} "${input_file}" >"${head_output}"
tail --lines ${tail_lines} "${head_output}" >"${output_file}"
rm -rf "${workdir}"

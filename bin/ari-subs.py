#!/usr/bin/env python
import re
import time
import os
import sys
import json
import argparse


def lj(fn):
    with open(fn) as handle:
        for line in handle:
            yield json.loads(line)


def gf(d, fn):
    return os.path.join(d, fn)


def timefmt(t, fmt):
    seconds = t % (24 * 3600)
    hours = int(seconds // 3600)
    seconds %= 3600
    minutes = int(seconds // 60)
    seconds %= 60
    (seconds, ms) = divmod(seconds, 1)
    seconds = int(seconds)
    ms = int(1000 * ms)

    if fmt == "webvtt":
        return f"{hours:02d}:{minutes:02d}:{seconds:02d}.{ms:03d}"
    else:
        return f"{hours:02d}:{minutes:02d}:{seconds:02d},{ms:03d}"


def script2timings(d):
    # Read in the sounds script
    with open(gf(d, "sounds.txt"), "r") as handle:
        script = handle.read()

    # This file contains something like the following:
    #
    #   file 'a46e0c84ae99c51a342c5b0dd51032ea.mp3'
    #   outpoint 1.283667
    #   file '68d15956eb04ebab450401be6501291c.mp3'
    #   outpoint 5.149833
    #
    # We need to find the ordering of the mp3s we want to smush. We will know the
    # timing from the json file.
    script = script.strip()
    ordering = re.findall(r"file '(.*).mp3'", script)
    timings = list(map(float, re.findall(r"duration (.*)", script)))
    return ordering, timings


def main(d, fmt="webvtt"):
    (ordering, timings) = script2timings(d)

    if fmt == "webvtt":
        print("WEBVTT\n\n")

    offset = 0
    for idx, (mid, timing) in enumerate(zip(ordering, timings)):
        if mid == "silence":
            offset += timing
            continue

        sound_data = list(lj(gf(d, mid + ".json")))
        sentences = [x for x in sound_data if x["type"] == "sentence"]

        sentence_offset = 0
        for si, sentence in enumerate(sentences):
            start = sentence["time"] / 1000
            # If it isn't the last sentence, we end when the next one starts.
            if si < len(sentences) - 1:
                end = sentences[si + 1]["time"] / 1000
            else:
                # Otherwise, it's the end of this section.
                end = timing

            print(idx)
            print(f"{timefmt(offset + start, fmt)} --> {timefmt(offset + end, fmt)}")
            print(sentence["value"])
            print()

        offset += timing


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert a video directory into subtitles"
    )
    parser.add_argument("directory", help="Directory")
    parser.add_argument(
        "--format", choices=["srt", "webvtt"], help="Choice of format", default="srt"
    )
    args = parser.parse_args()
    main(args.directory, fmt=args.format)

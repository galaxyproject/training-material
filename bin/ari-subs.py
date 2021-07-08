#!/usr/bin/env python
import re
import time
import os
import sys
import argparse


def timefmt(t, fmt):
    """
    Format a timestamp in the correct way for subtitles. This is the only
    reason we need a programming language and can't do this in bash easily.
    """
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
    with open(os.path.join(d, "sounds.txt"), "r") as handle:
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


def split_sentence(sentence, timing):
    split_sentence = sentence.split(' ')
    spaces = sentence.count(" ")
    parts = 1 + int(spaces // 20)
    part_off = int(spaces / parts)

    o = 0
    tacc = 0
    for i in range(parts):
        if i == parts - 1:
            yield (' '.join(split_sentence[o:]), tacc, tacc + timing / parts)
        else:
            yield (' '.join(split_sentence[o:o + part_off]), tacc, tacc + timing / parts)
        o += part_off
        tacc += timing / parts


def main(d, fmt="webvtt"):
    (ordering, timings) = script2timings(d)

    if fmt == "webvtt":
        print("WEBVTT\n\n")

    offset = 0
    idx = 0
    for (mid, timing) in zip(ordering, timings):
        # If there's a sentence, not silence
        if mid != 'silence':
            # load the text for the subtitles
            sound_data = os.path.join(d, mid + "-subtitle.txt")
            with open(sound_data, 'r') as handle:
                sentence = handle.read().strip()

            for (sen_part, time_prev, time_next) in split_sentence(sentence, timing):
                # print(sen_part, time_prev, time_next)
                # And print it.
                print(idx)
                print(f"{timefmt(offset + time_prev, fmt)} --> {timefmt(offset + time_next, fmt)}")
                print(sen_part)
                print()
            # sys.stderr.write(f'{sentence.count(" ")} {sentence}\n')

        offset += timing
        idx += 1


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

#!/usr/bin/env python
import os
import re
import subprocess
import tempfile
import json
import yaml
import sys

GTN_HOME = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")

ANSIBLE_HOST_OVERRIDE = "gat-1.be.training.galaxyproject.eu"
GTN_URL = "https://training.galaxyproject.org/training-material"
GXY_URL = "https://gat-1.be.training.galaxyproject.eu"
GIT_GAT = '/home/ubuntu/galaxy/'

GTN_URL = "http://localhost:4002/training-material"
GIT_GAT = '/home/hxr/arbeit/galaxy/git-gat'

data = yaml.safe_load(open(sys.argv[1], "r").read())
meta = data["meta"]
steps = data["steps"]

import os

if os.path.exists(".step.cache.json"):
    with open(".step.cache.json", "r") as handle:
        steps = json.load(handle)

# First pass, we'll record all of the audio, and nothing else.
for step in steps:
    # Ignore cached data
    if "mediameta" in step:
        continue

    # Need to record
    if "text" in step:
        tmp = tempfile.NamedTemporaryFile(mode="w", delete=False)
        tmp.write(step["text"])
        tmp.close()

        cmd = [
            "ruby",
            "bin/ari-synthesize.rb",
            "--aws",
            "-f",
            tmp.name,
            "--voice=" + meta["voice"]["id"],
            "--lang=" + meta["voice"]["lang"],
        ]
        if not meta["voice"]["neural"]:
            cmd.append("--non-neural")

        # Create the mp3
        mp3 = subprocess.check_output(cmd).decode("utf-8").strip()
        # And remove our temp file
        try:
            os.unlink(tmp)
        except:
            pass

        print(mp3)
        with open(mp3.replace(".mp3", ".json"), "r") as handle:
            mediameta = json.load(handle)
            del mediameta["value"]
        step["mediameta"] = mediameta
        step["mediameta"]["fn"] = mp3

    if "visual" not in step["data"]:
        raise Exception("Hey, missing a visual, we can't process this script.")

    with open(".step.cache.json", "w") as handle:
        json.dump(steps, handle)


def loadGitGatCommits(tutorial="admin/cvmfs"):
    subprocess.check_call(["git", "stash"], cwd=GIT_GAT)
    subprocess.check_call(["git", "checkout", "main"], cwd=GIT_GAT)
    results = subprocess.check_output(["git", "log", "--pretty=reference"], cwd=GIT_GAT).decode("utf-8")
    commits = results.strip().split("\n")[::-1]
    commits = [x for x in commits if tutorial in x]
    commitMap = {}
    for commit in commits:
        m = re.match("([0-9a-f]+) \(.*: (.*), [0-9-]*\)", commit)
        g = m.groups()
        commitMap[g[1]] = g[0]
    return commitMap


commitMap = loadGitGatCommits()


def runningGroup(steps):
    c = None
    for step in steps:
        if c is None:
            c = [step]
            continue

        # If the new visual is different, then
        if step["data"]["visual"] != c[-1]["data"]["visual"]:
            # Yield the current list
            yield c
            # And reset it with current step
            c = [step]
        else:
            c.append(step)
    yield c


def calculateSyncing(syncpoints, audio):
    """
    # This tracks where in the video precisely (mostly) the audio should start
    syncpoints = [
        {"msg": "checkout", "time": 162.959351},
        {"msg": "cmd", "time": 18899.63733},
        {"msg": "cmd", "time": 22978.593038},
    ]

    # Here are the lengths of the audio
    audio = [
        {"end": 13.608},
        {"end": 3.936},
        {"end": 7.488},
    ]

    So we need to get back something that has like (the commnd took 18.89
    seconds to run, but we have 13 seconds of audio, the next clip needs a 5
    second delay)

    {'end': 13.608} 162.959
    {'end': 3.936} 5128.677
    {'end': 7.488}  142.955
    """
    if len(syncpoints) != len(audio):
        print("????? Something odd is going on!")

    since = 0
    for (aud, syn) in zip(audio, syncpoints):
        delay = syn["time"] - since
        yield (int(delay), aud)
        since = syn["time"] + (aud["mediameta"]["end"] * 1000)


def muxAudioVideo(group, videoin, videoout, syncpoints):
    # We'll construct a complex ffmpeg filter here.
    #
    # There is incorrect quotation/line breaking in the filter_complex for clarity
    # ffmpeg -i video.mp4 \
    #   -i test3.mp3 \
    #   -i test3.mp3 \
    #   -filter_complex \
    #     "[1:a]atrim=start=0,adelay=500,asetpts=PTS-STARTPTS[a1];
    #      [2:a]atrim=start=0,adelay=20,asetpts=PTS-STARTPTS[a2];
    #      [a1][a2]concat=n=2:v=0:a=1[a]" \
    # -map 0:v -map "[a]" -y output.mp4
    #
    # We start by inputting all media (video, N audio samples)
    # Then we run a filter complex.
    # For each audio sample:
    #    we need to list it in the complex filter. We number then [1:a]... and refer to them as [a1]...
    #      atrim=start=0 to use the entire audio sample every time
    #                    http://ffmpeg.org/ffmpeg-filters.html#atrim
    #      adelay=500    value in miliseconds, we want to offset the first
    #                    clip relative to video start, the rest are offset
    #                    relative to each other (which in practice is a
    #                    negligible amount.)
    #                    http://ffmpeg.org/ffmpeg-filters.html#adelay
    #      asetpts=PTS-STARTPTS[..]
    #                    This is some magic I copied from SO. I don't super
    #                    get what it does, something about fixing
    #                    timestamps so there's no jitter in audio.
    #                    http://ffmpeg.org/ffmpeg-filters.html#asetpts
    # Then we list all audios with a concat filter and the number of samples.
    # [a1][a2]concat=n=2:v=0:a=1[a]
    # I read this as [audio1][audio2] concatentate filter with n=2samples : v(ideo)=off, a(udio)=on from [a] array of audio.
    # it works?
    #  and lastly we map things together into our output file.
    mux_cmd = ["ffmpeg", "-i", videoin]
    for g in group:
        mux_cmd.extend(["-i", g["mediameta"]["fn"]])
    mux_cmd.append("-filter_complex")

    delayResults = list(calculateSyncing(syncpoints, group))

    # The first audio sample must have a correct adelay for when that action happens.
    complex_filter = []
    concat_filter = ""
    for i2, (delay, g) in enumerate(delayResults, start=1):
        filt = f"[{i2}:a]atrim=start=0,adelay={delay},asetpts=PTS-STARTPTS[a{i2}]"
        print(filt)
        complex_filter.append(filt)
        concat_filter += f"[a{i2}]"

    final_concat = f"{concat_filter}concat=n={len(group)}:v=0:a=1[a]"
    complex_filter.append(final_concat)
    final_complex = ";".join(complex_filter)
    mux_cmd.extend([final_complex, "-map", "0:v", "-map", "[a]", videoout])
    # print(" ".join(mux_cmd))
    subprocess.check_call(mux_cmd)


def recordBrowser(idx):
    # Record our video. It'll start with a blank screen and page loading which we'll need to crop.
    if os.path.exists(f"video-{idx}.mp4"):
        return

    silent_video = f"video-{idx}-silent.mp4"
    silent_video_cropped = f"video-{idx}-cropped.mp4"
    cmd = [
        "/srv/galaxy/venv/bin/node",
        os.path.join(GTN_HOME, "bin", "video-browser-recorder.js"),
        f"scene-{idx}.json",
        silent_video,
    ]
    print(" ".join(cmd))
    resulting_script = json.loads(subprocess.check_output(cmd).decode("utf-8"))

    # Get the amount of time before the first scrollTo
    adelay = resulting_script[0]["time"]

    # Crop the 'init' portion of the video.
    cmd = [
        "ffmpeg",
        "-ss",
        f"{adelay}ms",
        "-i",
        silent_video,
        "-c",
        "copy",
        silent_video_cropped,
    ]
    print(" ".join(cmd))
    subprocess.check_call(cmd)

    # Build the video with sound.
    muxAudioVideo(group, silent_video_cropped, f"video-{idx}.mp4", resulting_script[1:])


def recordGtn(idx, group):
    # We've got N bits of text
    actions = [{"action": "goto", "target": GTN_URL + "/topics/admin/tutorials/cvmfs/tutorial.html"}]
    for g in group:
        actions.append(
            {"action": "scrollTo", "target": g["data"]["target"], "sleep": g["mediameta"]["end"],}
        )

    with open(f"scene-{idx}.json", "w") as handle:
        json.dump(actions, handle)

    recordBrowser(idx)


def recordGxy(idx, group):
    actions = [{"action": "goto", "target": GXY_URL,}]
    for g in group:
        action = {
            "action": g["data"]["action"],
            "target": g["data"]["target"],
            "value": g["data"].get("value", None),
            "sleep": g["mediameta"]["end"],
        }
        if action["action"] == "goto":
            action["target"] = GXY_URL + action["target"]

        actions.append(action)

    with open(f"scene-{idx}.json", "w") as handle:
        json.dump(actions, handle)

    recordBrowser(idx)


def recordTerm(idx, group):
    if os.path.exists(f"video-{idx}.mp4"):
        return

    actions = []
    for g in group:
        if "commit" in g["data"]:
            g["data"]["commitId"] = commitMap[g["data"]["commit"]]
            del g["code"]

        # t = g.get('mediameta', {'end': -1})['end']
        t = g["mediameta"]["end"]
        if "commitId" in g["data"]:
            actions.append({"action": "checkout", "time": t, "data": g["data"]["commitId"]})
        else:
            if "cmd" in g:
                cmd = g["cmd"]
            elif "cmd" in g["data"]:
                cmd = g["data"]["cmd"]
            else:
                print("????? SOMETHING IS WRONG")
            actions.append({"action": "cmd", "time": t, "data": cmd})
    with open(f"scene-{idx}.json", "w") as handle:
        json.dump(actions, handle)

    # Remove any previous versions of the cast.
    cast_file = f"{GTN_HOME}/scene-{idx}.cast"
    if os.path.exists(cast_file):
        os.unlink(cast_file)

    # Do the recording
    innercmd = [
        "bash",
        os.path.join(GTN_HOME, "bin", "video-term-recorder.sh"),
        f"{GTN_HOME}/scene-{idx}.json",
        f"{GTN_HOME}/scene-{idx}.log",
        ANSIBLE_HOST_OVERRIDE,
        GIT_GAT,
    ]
    cmd = [
        "asciinema",
        "rec",
        cast_file,
        "-t",
        f"Scene {idx}",
        "-c",
        " ".join(innercmd),
    ]
    print(' '.join(cmd))
    subprocess.check_call(cmd)
    # Convert to MP4
    subprocess.check_call(
        ["python", "asciicast2movie/asciicast2movie.py", f"scene-{idx}.cast", f"scene-{idx}.mp4",]
    )

    resulting_script = []
    with open(f"scene-{idx}.log", "r") as handle:
        for line in handle.readlines():
            line = line.strip().split("\t")
            resulting_script.append(
                {"time": float(line[0]) * 1000, "msg": line[1],}
            )

    # Mux with audio
    muxAudioVideo(group, f"scene-{idx}.mp4", f"video-{idx}.mp4", resulting_script)


# Next pass, we'll aggregate things of the same 'type'. This will make
# recording videos easier because we can more smoothly tween between steps.
# E.g. scrolling in GTN + waiting.  Or recording N things in the terminal and
# the audios for those.
for idx, group in enumerate(runningGroup(steps)):
    typ = group[0]["data"]["visual"]
    # print(typ, len(group), idx)
    if typ == "gtn":
        recordGtn(idx, group)
    elif typ == "terminal":
        recordTerm(idx, group)
    elif typ == "galaxy":
        recordGxy(idx, group)

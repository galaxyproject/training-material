#!/usr/bin/env python
import argparse
import glob
from collections import defaultdict
import os
import re
import shutil
import subprocess
import time
import yaml


def discover_trainings(topics_dir):
    """Auto-discover all topic metadata files."""
    for training in glob.glob(os.path.join(topics_dir, '*', 'metadata.yaml')):
        with open(training, 'r') as handle:
            training_data = yaml.load(handle)
            for material in training_data['material']:
                yield training.split('/')[-2], material['name'], material['title']


def safe_name(server, dashes=True):
    """Make human strings 'safe' for usage in paths."""
    safe_name = re.sub('\s', '_', server)
    if dashes:
        safe_name = re.sub('[^A-Za-z0-9_-]', '_', safe_name)
    else:
        safe_name = re.sub('[^A-Za-z0-9_]', '_', safe_name)

    return server


def get_badge_path(label, value, color):
    """Return a string representing the expected badge filename. Returns something like 'Training Name|Supported' or 'Training Name|Unsupported'."""
    safe_label = label.replace('@', '%40').replace(' ', '%20').replace('-', '--')
    safe_value = value.replace('@', '%40').replace(' ', '%20').replace('-', '--')
    return '%s-%s-%s.svg' % (safe_label, safe_value, color)


def realise_badge(badge, badge_cache_dir):
    """Download the badge to the badge_cache_dir (if needed) and return this real path to the user."""
    if not os.path.exists(os.path.join(badge_cache_dir, badge)):
        # Download the missing image
        cmd = [
            'wget', 'https://img.shields.io/badge/%s' % badge,
            '--quiet', '-O', os.path.join(badge_cache_dir, badge)
        ]
        subprocess.check_call(cmd)
        # Be nice to their servers
        time.sleep(1)
    return os.path.join(badge_cache_dir, badge)


def badge_it(label, value, color, CACHE_DIR, instance, identifier, output_dir):
    # Get a path to a (cached) badge file.
    real_badge_path = realise_badge(get_badge_path(
        label, value, color
    ), CACHE_DIR)
    # Deteremine the per-instance output name
    output_filename = safe_name(instance) + '__' + safe_name(identifier, dashes=True) + '.svg'
    output_filepath = os.path.join(args.output, output_filename)
    # Copy the badge to a per-instance named .svg file.
    shutil.copy(real_badge_path, output_filepath)
    return output_filename


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build the badge directory for instances to use.')
    parser.add_argument('--public-server-list', help='Url to access the public galaxy server list at',
                        default='https://raw.githubusercontent.com/martenson/public-galaxy-servers/master/servers.csv')
    parser.add_argument('--topics-directory', help='Path to the topics directory', default='./topics/')
    parser.add_argument('--instances', help='File containing the instances and their supported trainings', default='metadata/instances.yaml')

    parser.add_argument('--output', help='Path to the the directory where the badges should be stored. The directory will be created if it does not exist.', default='output')
    args = parser.parse_args()

    # Validate training dir argument
    if not os.path.exists(args.topics_directory) and os.path.is_dir(args.topics_directory):
        raise Exception("Invalid topics directory")
    all_trainings = list(discover_trainings(args.topics_directory))
    trainings = {x[1]: x[2] for x in all_trainings}

    topic_counts = defaultdict(int)
    for (topic, _, _) in all_trainings:
        topic_counts[topic] += 1

    training_keys = sorted(trainings.keys())
    if len(trainings) == 0:
        raise Exception("No trainings discovered!")

    # Create output directory if not existing.
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Also check/create the badge cache directory.
    CACHE_DIR = os.path.join(args.output, '.cache')
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)

    # Load the validated list of instances which support trainings
    with open(args.instances, 'r') as handle:
        data = yaml.load(handle)

    # Collect a list of instances seen
    instances = []
    for topic in data:
        for training in data[topic]:
            for instance in data[topic][training]:
                data[topic][training][instance]['supported'] = True
                instances.append(instance)
    instances = sorted(set(instances))

    # Mark the unsupported ones as such for easier processing later.
    for topic in data:
        for training in data[topic]:
            for instance in instances:
                # Not in one of the existing supported ones
                if instance not in data[topic][training]:
                    data[topic][training][instance] = {'supported': False}

    index_html = open(os.path.join(args.output, 'index.html'), 'w')

    index_html.write("""
<html><head></head>
<body>
<h1>Galaxy Training Support</h1>
<p>
These badges are the results of regularly checking in with your Galaxy
instance to see if you support the various training courses that are available
to end users. This requires having all of the appropriate tools installed and
possibly datasets in specifically named data libraries.
</p>""")

    # Map of instance -> badges
    instance_badges = {}
    # Count of tutorials in each topic.
    for topic in data:
        # All trainings, not just those available
        for training in sorted(data[topic]):
            for instance in data[topic][training]:
                if instance not in instance_badges:
                    instance_badges[instance] = {}

                if topic not in instance_badges[instance]:
                    instance_badges[instance][topic] = []

                # If available, green badge
                is_supported = data[topic][training][instance]['supported']
                # We'll only place the badge in the HTML if the training is
                # supported (but the unavailable badge will still be available
                # in case they ever go out of compliance.)
                if is_supported:
                    output_filename = badge_it(
                        trainings[training],
                        'Supported', 'green',
                        CACHE_DIR, instance, training, args.output
                    )
                    instance_badges[instance][topic].append(output_filename)
                else:
                    badge_it(
                        trainings[training],
                        'Unsupported', 'lightgrey',
                        CACHE_DIR, instance, training, args.output
                    )


    # All instances, not just checked
    for instance in sorted(instance_badges):
        total = sum([len(instance_badges[instance][topic]) for topic in instance_badges[instance]])

        if total == 0:
            continue

        index_html.write('<h2 id="' + safe_name(instance, dashes=True) + '">' + instance + '</h2>')
        index_html.write('<h3>Per-Training Badge</h3>')
        index_html.write('<ul>')
        for topic in instance_badges[instance]:
            for badge in instance_badges[instance][topic]:
                index_html.write('<li><img src="' + badge + '"/></li>')
        index_html.write('</ul>')
        index_html.write('<h3>Training Group Badges</h3>')
        index_html.write('<ul>')
        for topic in instance_badges[instance]:
            # Get the number of badges in this topic.
            count = len(instance_badges[instance][topic])

            if float(count) / topic_counts[topic] > 0.90:
                color = 'green'
            elif float(count) / topic_counts[topic] > 0.25:
                color = 'orange'
            else:
                color = 'red'

            output_filename = badge_it(
                topic, '%s%%2f%s' % (count, topic_counts[topic]), color,
                CACHE_DIR, instance, topic, args.output
            )

            if count > 0:
                index_html.write('<li><img src="' + output_filename + '"/></li>')
        index_html.write('</ul>')

    index_html.write("</body></html>")

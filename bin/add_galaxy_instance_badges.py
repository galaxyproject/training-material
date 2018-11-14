#!/usr/bin/env python
import argparse
import glob
from collections import defaultdict
import os
import re
import subprocess
import time
import yaml
DRY_RUN = False


def discover_trainings(topics_dir):
    """Auto-discover all topic metadata files."""
    for training_dir in glob.glob(os.path.join(topics_dir, '*')):

        with open(os.path.join(training_dir, 'metadata.yaml'), 'r') as handle:
            training_data = yaml.load(handle)

        training = {
            'title': training_data['title'],
            'trainings': {},
        }

        for material in glob.glob(os.path.join(training_dir, 'tutorials', '*', 'tutorial.md')) + glob.glob(os.path.join(training_dir, 'tutorials', '*', 'slides.html')):
            with open(material, 'r') as handle:
                material_data = yaml.load_all(handle)
                material_data = next(material_data)

            name = material.split('/')[-2]
            training['trainings'][name] = material_data['title']

        training['count'] = len(training['trainings'].keys())

        yield training_data['name'], training


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
    print(label, value, color)
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
        if not DRY_RUN:
            subprocess.check_call(cmd)
            time.sleep(1)
        else:
            print(' '.join(cmd))
        # Be nice to their servers
    return os.path.join(badge_cache_dir, badge)


def badge_it(label, value, color, CACHE_DIR, identifier_parts, output_dir):
    # Get a path to a (cached) badge file.
    real_badge_path = realise_badge(get_badge_path(
        label, value, color
    ), CACHE_DIR)
    # Deteremine the per-instance output name
    output_filedir = os.path.join(args.output, *map(safe_name, identifier_parts[0:-1]))
    if not os.path.exists(output_filedir):
        os.makedirs(output_filedir)

    output_filename = safe_name(identifier_parts[-1]) + '.svg'
    # Ensure dir exists
    output_filepath = os.path.join(output_filedir, output_filename)

    # Copy the badge to a per-instance named .svg file.
    up = ['..'] * (len(identifier_parts) - 1)
    symlink_source = os.path.join(*up, real_badge_path[len('badges/'):])
    if not DRY_RUN:
        # Remove it if it exists, since this is easier than testing for
        # equality.
        if os.path.exists(output_filepath):
            os.unlink(output_filepath)

        # Now (re-)create the symlink
        os.symlink(symlink_source, output_filepath)
    else:
        print(' '.join(['ln -s ', symlink_source, output_filepath]))
    return output_filename


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build the badge directory for instances to use.')
    parser.add_argument('--public-server-list', help='Url to access the public galaxy server list at',
                        default='https://raw.githubusercontent.com/martenson/public-galaxy-servers/master/servers.csv')
    parser.add_argument('--topics-directory', help='Path to the topics directory', default='./topics/')
    parser.add_argument('--instances', help='File containing the instances and their supported trainings', default='metadata/instances.yaml')

    parser.add_argument('--output', help='Path to the the directory where the badges should be stored. The directory will be created if it does not exist.', default='badges')
    args = parser.parse_args()

    # Validate training dir argument
    if not os.path.exists(args.topics_directory) and os.path.is_dir(args.topics_directory):
        raise Exception("Invalid topics directory")
    all_trainings = {k: v for (k, v) in discover_trainings(args.topics_directory)}

    # Create output directory if not existing.
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Also check/create the badge cache directory.
    CACHE_DIR = os.path.join(args.output, 'cache')
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)

    # Load the validated list of instances which support trainings
    with open(args.instances, 'r') as handle:
        data = yaml.load(handle)

    # Collect a list of instances seen
    instances = []
    for topic in data:
        for training in data[topic]['tutorials']:
            for instance in data[topic]['tutorials'][training]['instances']:
                data[topic]['tutorials'][training]['instances'][instance]['supported'] = True
                instances.append(instance)
    # All of these instances support at least one training.
    instances = sorted(set(instances))

    # Mark the unsupported ones as such for easier processing later.
    for topic in data:
        for training in data[topic]['tutorials']:
            for instance in instances:
                # Not in one of the existing supported ones
                if instance not in data[topic]['tutorials'][training]['instances']:
                    data[topic]['tutorials'][training]['instances'][instance]['supported'] = False

    # Map of instance -> badges
    instance_badges = {}

    # Count of tutorials in each topic.
    for topic in data:
        # All trainings, not just those available
        for training in sorted(data[topic]['tutorials']):
            for instance in data[topic]['tutorials'][training]['instances']:
                if instance not in instance_badges:
                    instance_badges[instance] = {}

                if topic not in instance_badges[instance]:
                    instance_badges[instance][topic] = []

                # If available, green badge
                is_supported = data[topic]['tutorials'][training]['instances'][instance]['supported']

                # We'll only place the badge in the HTML if the training is
                # supported (but the unavailable badge will still be available
                # in case they ever go out of compliance.)

                label = all_trainings[topic]['trainings'][training]
                if is_supported:
                    output_filename = badge_it(
                        label,
                        'Supported', 'green',
                        CACHE_DIR, (instance, topic, training), args.output
                    )
                    instance_badges[instance][topic].append(output_filename)
                else:
                    badge_it(
                        label,
                        'Unsupported', 'lightgrey',
                        CACHE_DIR, (instance, topic, training), args.output
                    )

    # All instances, not just checked
    for instance in sorted(instance_badges):
        total = sum([len(instance_badges[instance][topic]) for topic in instance_badges[instance]])

        if total == 0:
            continue

        for topic in instance_badges[instance]:
            # Get the number of badges in this topic.
            count = len(instance_badges[instance][topic])

            if float(count) / all_trainings[topic]['count'] > 0.90:
                color = 'green'
            elif float(count) / all_trainings[topic]['count'] > 0.25:
                color = 'orange'
            else:
                color = 'red'

            output_filename = badge_it(
                all_trainings[topic]['title'], '%s%%2f%s' % (count, all_trainings[topic]['count']), color,
                CACHE_DIR, (instance, topic), args.output
            )

#!/usr/bin/env python

import os
import argparse
import yaml
from bioblend import galaxy
import pprint

log_dir = "logs"

def check_instance(galaxy_url, tool_filepath, log_file):
    """ Check if a given tutorial can be run on a given Galaxy instance """
    gi = galaxy.GalaxyInstance(url=galaxy_url)
    # find tools needed for the tutorial
    with open(tool_filepath) as f:
        tool_yaml = yaml.safe_load(f)
    tools_required = {}
    for tool in tool_yaml['tools']:
        tools_required.setdefault(
            tool['name'],
            {'found': False, 'revisions': {}})
        if 'revisions' in tool:
            for revision in tool['revisions']:
                tools_required[tool['name']]['revisions'].setdefault(revision, False)
    # check if all tools are installed
    tools_list = gi.tools.get_tools()
    for t in tools_list:
        if 'tool_shed_repository' not in t:
            continue
        name = t['tool_shed_repository']['name']
        revision = t['tool_shed_repository']['changeset_revision']
        if name in tools_required:
            tools_required[name]['found'] = True
            if revision in tools_required[name]['revisions']:
                tools_required[name]['revisions'][revision] = True
    # write in log file what is installed and what is missing
    tool_requirements_satisfied = True
    with open(log_file, 'a') as f:
        f.write("Tutorial: %s\n" % (tool_filepath))
        for tool in tools_required:
            f.write(" - %s: %s\n" % (tool, tools_required[tool]['found']))
            tool_requirements_satisfied &= tools_required[tool]['found']
            for revision in tools_required[tool]['revisions']:
                f.write("  - %s: %s\n" % (revision, tools_required[tool]['revisions'][revision]))
        f.write("\n")
    return tool_requirements_satisfied


def check_tutorial(tool_file, public_galaxy_servers):
    """ Check which public Galaxy servers can run the tool in a tutorial """
    working_galaxy_servers = []
    for server in public_galaxy_servers:
        log_file = os.path.join(log_dir, server['name'])
        can_run = check_instance(server['url'], tool_file, log_file)
        if can_run:
            working_galaxy_servers.append(server)


def check_tutorials(public_galaxy_servers):
    """ Extract a dictionary with for each topic the list tutorials to tests """
    to_test = {}
    os.makedirs(log_dir, exist_ok=True)
    # remove existing log files
    for f in os.listdir(log_dir):
        os.remove(os.path.join(log_dir, f))
    # parse topics
    for topic in os.listdir("topics"):
        to_test.setdefault(topic, [])
        topic_path = os.path.join("topics", topic)
        tutorials_path = os.path.join(topic_path, "tutorials")
        # load topic metadata
        metadata_file = os.path.join(topic_path, "metadata.yaml")
        with open(metadata_file, 'r') as f:
            topic_metadata = yaml.safe_load(f)
        print(topic_metadata)
        # parse the tutorials of a topic
        for tutorial in os.listdir(tutorials_path):
            tutorial_path = os.path.join(tutorials_path, tutorial)
            # pass if the tools.yaml file does not exist
            tool_file = os.path.join(tutorial_path, 'tools.yaml')
            if not os.path.exists(tool_file):
                continue
            # find the Galaxy servers with the needed tools
            #working_galaxy_servers = check_tutorial(tool_file, public_galaxy_servers)
            # add the information to the metadata of the tutorial
        # write the metadata
        pp = pprint.PrettyPrinter(indent=4)
        with open(metadata_file, 'w') as f:
            yaml.safe_dump(topic_metadata, f)
        #pp.pprint(topic_metadata)
    return to_test


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract which public Galaxy servers can run which tutorial and add the information in the metadata')
    parser.add_argument(
        '-p',
        '--public_galaxy_servers',
        default="galaxy_instances.yaml",
        help="Path to the YAML file with the list of public Galaxy servers")
    args = parser.parse_args()

    with open(args.public_galaxy_servers) as f:
        public_galaxy_servers = yaml.safe_load(f)

    check_tutorials(public_galaxy_servers)

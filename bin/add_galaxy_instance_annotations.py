#!/usr/bin/env python

import os
import argparse
import yaml
from bioblend import galaxy
import requests
import io
import pandas as pd

log_dir = "logs"


def extract_public_galaxy_servers():
    """ Extract a pandas Data frame with the public Galaxy servers """
    r = requests.get('https://raw.githubusercontent.com/martenson/public-galaxy-servers/master/servers.csv')
    f = io.StringIO(r.text)
    public_galaxy_servers = pd.read_csv(f, sep=",")
    return public_galaxy_servers


def test_server_connection(galaxy_url):
    """ Test if the server is accessible and return the GalaxyInstance object """
    gi = galaxy.GalaxyInstance(url=galaxy_url)
    try:
        gi.config.get_config()
    except:
        print("Can not connect to %s" % galaxy_url)
        return None
    return gi


def check_instance(gi, tool_filepath, log_file):
    """ Check if a given tutorial can be run on a given Galaxy instance """
    # get tools in the instance
    try:
        tools_list = gi.tools.get_tools()
    except:
        print("Can not geet tools from %s" % galaxy_url)
        return None
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
    to_conserve = []
    to_drop = []
    working_galaxy_servers = None
    # parse the public servers
    for index, server in public_galaxy_servers.iterrows():
        # test if server is accessible and remove from the list if not
        gi = test_server_connection(server['url'])
        if gi is None:
            to_drop.append(index)
            continue
        # check if the tools are accessible on the serve
        log_file = os.path.join(log_dir, server['name'])
        can_run = check_instance(gi, tool_file, log_file)
        if can_run is None:
            to_drop.append(index)
        elif can_run:
            to_conserve.append(index)
    # extract the info of public server that can be used for the tutorial
    if len(to_conserve) > 0:
        working_galaxy_servers = public_galaxy_servers.ix[to_conserve]
    # remove the non accessible servers
    public_galaxy_servers = public_galaxy_servers.drop(to_drop)
    return working_galaxy_servers, public_galaxy_servers


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
        #print(topic_metadata)
        # parse the tutorials of a topic
        for tutorial in os.listdir(tutorials_path):
            print(tutorial)
            tutorial_path = os.path.join(tutorials_path, tutorial)
            # pass if the tools.yaml file does not exist
            tool_file = os.path.join(tutorial_path, 'tools.yaml')
            if not os.path.exists(tool_file):
                continue
            # find the Galaxy servers with the needed tools
            working_galaxy_servers, public_galaxy_servers = check_tutorial(tool_file, public_galaxy_servers)
            print(working_galaxy_servers)
            # add the information to the metadata of the tutorial
        # write the metadata
        with open(metadata_file, 'w') as f:
            yaml.safe_dump(topic_metadata, f)
        #pp.pprint(topic_metadata)
    return to_test


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract which public Galaxy servers can run which tutorial and add the information in the metadata')
    args = parser.parse_args()

    public_galaxy_servers = extract_public_galaxy_servers()
    check_tutorials(public_galaxy_servers)

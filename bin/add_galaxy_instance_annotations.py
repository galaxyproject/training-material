#!/usr/bin/env python

import os
import argparse
import yaml
from bioblend import galaxy


def check_instance(galaxy_url, tutorial_path):
    """ Check if a given tutorial can be run on a given Galaxy instance """
    # TODO: support specific revision checking
    # TODO: discussion: for suites, encourage to list all tools actually needed in yaml, not just the suite?

    gi = galaxy.GalaxyInstance(url=galaxy_url)

    # find tools needed for the tutorial
    with open(os.path.join(tutorial_path, 'tools.yaml')) as f:
        tool_yaml = yaml.safe_load(f)

    tools_required = []
    for tool in tool_yaml['tools']:
        tools_required.append(tool['name'])

    # check if all tools are installed
    tools_list = gi.tools.get_tools()
    tools_available = []
    for t in tools_list:
        try:
            tools_available.append(t['tool_shed_repository']['name'])
        except KeyError:
            pass

    tool_requirements_satisfied = set(tools_required) < set(tools_available)

    # output what is installed and what is missing
    print("This Galaxy instance ("+ galaxy_url + ") can run my tutorial: " + str(tool_requirements_satisfied))
    for tool in tools_required:
        print(" - " + tool + ": " + ("present" if tool in tools_available else "not present"))

    return tool_requirements_satisfied


def extract_tutorials_to_test():
    """ Extract a dictionary with each topic and tutorials for each topics """
    to_test = {}
    for topic in os.listdir("topics"):
        print(topic)
        to_test.setdefault(topic, [])
        topic_path = os.path.join("topics", topic)
        tutorials_path = os.path.join(topic_path, "tutorials")
        for tutorial in os.listdir(tutorials_path):
            print("  %s" % tutorial)
            tutorial_path = os.path.join(tutorials_path, tutorial)
            to_test[topic].append(tutorial)
    print(to_test)


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
    print(public_galaxy_servers)

    tutorial_list = extract_tutorials_to_test()
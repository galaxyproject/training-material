#!/usr/bin/env python

import os
import argparse
import yaml
from bioblend import galaxy


def check_instance(gi, tutorial_path):
    """ Check if a given tutorial can be run on a given Galaxy instance """
    # TODO: support specific revision checking

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
    print "This galaxy instance can run my tutorial: " + str(tool_requirements_satisfied)

    for tool in tools_required:
        print " - " + tool + ": " + ("present" if tool in tools_available else "not present")

    return tool_requirements_satisfied



def check_from_metadata(gi, tutorial):
    """ Check the instances mentioned in the metadata file of a given tutorial """
    pass

def check_repo(gi):
    """ Check instances in metadata files of all tutorials in the repo """
    pass



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check if tools required for a tutorial are present on a given Galaxy instance.')
    parser.add_argument('-c', '--command', choices=['check_from_url', 'check_from_metadata', 'check_repo'], required=True)
    parser.add_argument('-g', '--galaxy_url', required=True)
    parser.add_argument('-t', '--tutorial', help="Path to tutorial directory")
    parser.add_argument('-m', '--mapping_yaml', default="galaxy_instances.yaml",
                        help="location of yaml file with mappings of Galaxy instance names to urls")

    args = parser.parse_args()

    gi = galaxy.GalaxyInstance(url=args.galaxy_url)

    if args.command == 'check_from_url':
        check_instance(gi, args.tutorial)
    elif args.command == 'check_from_metadata':
        check_from_metadata(gi, args.tutorial)
    elif args.command == 'check_repo':
        check_repo(gi)

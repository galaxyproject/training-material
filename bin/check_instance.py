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
    print "This Galaxy instance ("+ galaxy_url + ") can run my tutorial: " + str(tool_requirements_satisfied)

    for tool in tools_required:
        print " - " + tool + ": " + ("present" if tool in tools_available else "not present")

    return tool_requirements_satisfied



def check_from_metadata(tutorial_path, mapping_yaml_file):
    """ Check the instances mentioned in the metadata file of a given tutorial """

    '''
    idea is for each tutorial to have a list of galaxy instances able to run the tutorial in the metadata

    in tutorial's metadata.yaml:
        galaxy_instances:
            - name: usegalaxy
            - name: freiburg

    in galaxy_instances.yaml in repo root:

    galaxies:
        - name: usegalaxy
          url: https://usegalaxy.org/
        - name: freiburg
          url: http://galaxy.uni-freiburg.de/
    '''

    # load tutorial metadata file
    with open(os.path.join(tutorial_path, 'metadata.yaml')) as f:
        metadata_yaml = yaml.safe_load(f)

    # load galaxy instances mapping file
    with open(mapping_yaml_file) as f:
        mapping_yaml = yaml.safe_load(f)

    # check if tools required for tutorial are present on each of the listed Galaxy instances
    for instance in metadata_yaml['galaxy_instances']:
        instance_name = instance['name']
        instance_url = get_galaxy_url(mapping_yaml, instance_name)

        if instance_url != None:
            check_instance(instance_url, tutorial_path)


def check_repo(gi):
    """ Check instances in metadata files of all tutorials in the repo """
    pass

def get_galaxy_url(mapping_yaml, name):
    for galaxy in mapping_yaml['galaxies']:
        if galaxy['name'] == name:
            return galaxy['url']

    print "WARNING: Galaxy instance '"+name+"' not found in mapping file, please correct"
    return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check if tools required for a tutorial are present on a given Galaxy instance.')
    parser.add_argument('-c', '--command', choices=['check_from_url', 'check_from_metadata', 'check_repo'], required=True)
    parser.add_argument('-g', '--galaxy_url')
    parser.add_argument('-t', '--tutorial', help="Path to tutorial directory")
    parser.add_argument('-m', '--mapping_yaml', default="galaxy_instances.yaml",
                        help="location of yaml file with mappings of Galaxy instance names to urls")

    args = parser.parse_args()

    if args.command == 'check_from_url':
        check_instance(args.galaxy_url, args.tutorial)
    elif args.command == 'check_from_metadata':
        check_from_metadata(args.tutorial, args.mapping_yaml)
    elif args.command == 'check_repo':
        check_repo(gi)

#!/usr/bin/env python
import argparse
import csv
import io
import json
import logging
import multiprocessing.pool
import os
import requests
import yaml


def extract_public_galaxy_servers():
    """Extract a pandas Data frame with the public Galaxy servers"""
    r = requests.get('https://raw.githubusercontent.com/martenson/public-galaxy-servers/master/servers.csv')
    f = io.StringIO(r.text)
    reader = csv.reader(f, delimiter=",")
    header = next(reader)
    for row in reader:
        yield dict(zip(header, row))


def fetch_and_extract_individual_server_tools(server):
    # request the tools via the API
    url = '%s/api/tools' % server['url'].rstrip('/')
    try:
        response = requests.get(url, timeout=20)
    except:
        print(server['name'] + " Error of connection")
        return

    # check status
    if response.status_code != requests.codes.ok:
        print(server['name'] + " Bad status (%s)" % response.status_code)
        return

    # check content
    if response.text.find("</html>") != -1:
        print(server['name'] + " No JSON output")
        return

    # extract the list of tools in this instance
    found_tools = set()
    try:
        response_json = response.json()
    except json.decoder.JSONDecodeError:
        print(server['name'] + " Invalid JSON")
        return

    for section in response_json:
        if 'elems' in section:
            for tool in section['elems']:
                found_tools.add('/'.join( tool['id'].split('/')[:4] ))

    return server['name'], {
        'url': server['url'],
        'tools': set(found_tools)
    }


def extract_public_galaxy_servers_tools():
    """Extract the tools from the public Galaxy servers using their API"""
    server_tools = {}

    to_process = []
    for server in extract_public_galaxy_servers():
        to_process.append(server)

    pool = multiprocessing.pool.ThreadPool(processes=20)
    processed = pool.map(fetch_and_extract_individual_server_tools, to_process, chunksize=1)
    pool.close()

    for server_data in processed:
        if server_data:
            server_tools[server_data[0]] = server_data[1]

    return server_tools


def check_tutorial(topic, tutorial, tool_filepath, server_tools):
    """Check which public Galaxy servers can run the tool in a tutorial"""

    # Extract the list of tools for the tutorial
    with open(tool_filepath, "r") as f:
        tool_yaml = yaml.safe_load(f)

    tutorial_tools = set()
    for tool in tool_yaml['tools']:
        tutorial_tools.add('toolshed.g2.bx.psu.edu/repos/%s/%s' % (tool['owner'], tool['name']))

    servers_supported = {}
    servers_unsupported = {}

    for server in server_tools:
        # check overlap
        if tutorial_tools.issubset(server_tools[server]['tools']):
            servers_supported[server] = server_tools[server]['url']
        else:
            missing_tools = tutorial_tools - set(server_tools[server]['tools'])
            logging.debug('MISSING:%s:%s/%s: %s', server, topic, tutorial, ', '.join(missing_tools))
            servers_unsupported[server] = server_tools[server]['url']

    return servers_supported, servers_unsupported


def check_tutorials(server=None):
    """Check all tutorials to extract the instances on which the tutorials can
    be run and add this information to metadata/instances.yaml file"""
    if server:
        name, server = fetch_and_extract_individual_server_tools({'url': server, 'name': 'none'})
        server_tools = {name: server}
    else:
        if os.path.exists('.cache.yaml'):
            with open('.cache.yaml', 'r') as handle:
                server_tools = yaml.load(handle)
        else:
            server_tools = extract_public_galaxy_servers_tools()
            with open('.cache.yaml', 'w') as handle:
                yaml.dump(server_tools, handle)

    instance_annot = {}
    for topic in os.listdir("topics"):
        topic_dir = os.path.join("topics", topic)
        tutos_dir = os.path.join(topic_dir, "tutorials")
        instance_annot.setdefault(topic, {'supported': False, 'tutorials': {}})
        if not os.path.exists(tutos_dir):
            continue

        for tutorial in os.listdir(tutos_dir):
            # extract the tool file
            tool_filepath = os.path.join(tutos_dir, tutorial, "tools.yaml")
            if not os.path.exists(tool_filepath):
                continue
            # extract the instances on which the tutorial can be run
            supported, unsupported = check_tutorial(topic, tutorial, tool_filepath, server_tools)

            # Training-focused view
            instance_annot[topic]['tutorials'].setdefault(tutorial, {'supported': False, 'instances': {}})
            instance_annot[topic]['tutorials'][tutorial]['instances'] = {k: {'url': v, 'supported': True} for (k, v) in supported.items()}
            instance_annot[topic]['tutorials'][tutorial]['supported'] = any(instance_annot[topic]['tutorials'][tutorial]['instances'])
            instance_annot[topic]['supported'] |= instance_annot[topic]['tutorials'][tutorial]['supported']
            instance_annot[topic]['tutorials'][tutorial]['instances'].update({k: {'url': v, 'supported': False} for (k, v) in unsupported.items()})

    # If we're checking a single server, do not update the metadata file.
    if server:
        return

    # add the conserved servers to the metadata/instances.yaml file
    inst_file = os.path.join("metadata", "instances.yaml")
    with open(inst_file, "w") as f:
        yaml.dump(instance_annot, f, default_flow_style=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract which public Galaxy servers can run for the tutorials and add this information to a instance file')
    parser.add_argument("--verbose", action='store_true')
    parser.add_argument("--server", help='Only check a single server, skip fetching of the public servers list.')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    check_tutorials(server=args.server)

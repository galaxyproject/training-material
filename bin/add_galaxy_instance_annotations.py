#!/usr/bin/env python
import json
import os
import argparse
import yaml
import requests
import io
import pandas as pd


def extract_public_galaxy_servers():
    """Extract a pandas Data frame with the public Galaxy servers"""
    r = requests.get('https://raw.githubusercontent.com/martenson/public-galaxy-servers/master/servers.csv')
    f = io.StringIO(r.text)
    public_galaxy_servers = pd.read_csv(f, sep=",")
    return public_galaxy_servers


def extract_public_galaxy_servers_tools():
    """Extract the tools from the public Galaxy servers using their API"""
    servers = extract_public_galaxy_servers()
    server_tools = {}
    for index, server in servers.iterrows():
        print(server['name'])
        # request the tools via the API
        url = '%s/api/tools' % server['url'].rstrip('/')
        try:
            response = requests.get(url, timeout=20)
        except:
            print("Error of connection")
        # check status
        if response.status_code != requests.codes.ok:
            print("Bad status (%s)" % response.status_code)
            continue
        # check content
        if response.text.find("</html>") != -1:
            print("No JSON output")
            continue
        # extract the list of tools in this instance
        found_tools = set()
        try:
            response_json = response.json()
        except json.decoder.JSONDecodeError:
            print("Invalid JSON")
            continue

        for section in response_json:
            if 'elems' in section:
                for tool in section['elems']:
                    found_tools.add('/'.join( tool['id'].split('/')[:4] ))
        # save the server with its tools
        server_tools[server['name']] = {
            'url': server['url'],
            'tools': found_tools
        }
    return server_tools


def extract_tool_list(tool_filepath):
    """Extract the list of tools for the tutorial"""
    with open(tool_filepath, "r") as f:
        tool_yaml = yaml.safe_load(f)
    tools = set()
    for tool in tool_yaml['tools']:
        tools.add('toolshed.g2.bx.psu.edu/repos/%s/%s' % (tool['owner'], tool['name']))
    return tools


def check_tutorial(tool_filepath, serv_tools):
    """Check which public Galaxy servers can run the tool in a tutorial"""
    exp_tools = extract_tool_list(tool_filepath)
    # parse the public servers
    cons_servers = {}
    for server in serv_tools:
        # check overlap
        res = exp_tools.issubset(serv_tools[server]['tools'])
        if res:
            cons_servers[server] = {'url': serv_tools[server]['url']}
    return cons_servers


def check_tutorials():
    """Check all tutorials to extract the instances on which the tutorials can
    be run and add this information to metadata/instances.yaml file"""
    server_tools = extract_public_galaxy_servers_tools()
    instance_annot = {}
    for topic in os.listdir("topics"):
        topic_dir = os.path.join("topics", topic)
        tutos_dir = os.path.join(topic_dir, "tutorials")
        for tutorial in os.listdir(tutos_dir):
            # extract the tool file
            tool_filepath = os.path.join(tutos_dir, tutorial, "tools.yaml")
            if not os.path.exists(tool_filepath):
                continue
            # extract the instances on which the tutorial can be run
            cons_servers = check_tutorial(tool_filepath, server_tools)
            # add the conserved servers to the main dictionary if not empty
            if not cons_servers:
                continue
            instance_annot.setdefault(topic, {})
            instance_annot[topic][tutorial] = cons_servers
    # add the conserved servers to the metadata/instances.yaml file
    inst_file = os.path.join("metadata", "instances.yaml")
    with open(inst_file, "w") as f:
        yaml.dump(instance_annot, f, default_flow_style=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract which public Galaxy servers can run for the tutorials and add this information to a instance file')
    args = parser.parse_args()
    check_tutorials()

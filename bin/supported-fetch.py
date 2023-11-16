#!/usr/bin/env python
import argparse
import csv
import io
import json
import logging
import multiprocessing.pool
import os
import requests


def fetch_and_extract_individual_server_tools(server):
    # request the tools via the API
    url = '%s/api/tools?in_panel=False' % server['url'].rstrip('/')
    try:
        response = requests.get(url, timeout=20)
    except:
        print(server['name'] + " Connection Timeout (20s)")
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
    try:
        response_json = response.json()
    except json.decoder.JSONDecodeError:
        print(server['name'] + " Invalid JSON")
        return

    found_tools = set()
    for tool in response_json:
        found_tools.add(tool['id'])

    return server['name'], {
        'url': server['url'],
        'tools': list(set(found_tools))
    }


def extract_public_galaxy_servers_tools():
    """Extract the tools from the public Galaxy servers using their API"""
    server_tools = {}

    to_process = []
    serverlist = requests.get('https://galaxyproject.org/use/feed.json').json()
    for server in serverlist:
        # We intentionally drop all usegalaxy.eu subdomains. They're all the
        # same as the top level domain and just pollute the supported instances
        # list.
        if '.usegalaxy.eu' in server['url']:
            continue
        # Apparently the french do it too
        if '.usegalaxy.fr' in server['url']:
            continue
        # The aussies will soon
        if '.usegalaxy.org.au' in server['url']:
            continue
        # No test servers permitted
        if 'test.' in server['url']:
            continue

        s = { 'name': server['title'], 'url': server['url'] }
        to_process.append(s)

    pool = multiprocessing.pool.ThreadPool(processes=20)
    processed = pool.map(fetch_and_extract_individual_server_tools, to_process, chunksize=1)
    pool.close()

    for server_data in processed:
        if server_data:
            server_tools[server_data[0]] = server_data[1]

    return server_tools


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract which public Galaxy servers can run specific tools')
    args = parser.parse_args()

    server_tools = extract_public_galaxy_servers_tools()
    # Reverse the mapping
    tool_servers = {
        'servers': [],
        'tools': {},
    }
    for idx, (server_name, server_data) in enumerate(server_tools.items()):
        tool_servers['servers'].append({
            'url': server_data['url'],
            'name': server_name,
        })
        for tool in server_data['tools']:
            if tool.count('/') > 4:
                tool_id = '/'.join(tool.split('/')[:5])
                tool_version = '/'.join(tool.split('/')[5:])
                if tool_id not in tool_servers['tools']:
                    tool_servers['tools'][tool_id] = {}
                if tool_version not in tool_servers['tools'][tool_id]:
                    tool_servers['tools'][tool_id][tool_version] = []

                tool_servers['tools'][tool_id][tool_version].append(idx)
            else:
                if tool not in tool_servers['tools']:
                    tool_servers['tools'][tool] = {"_": []}
                tool_servers['tools'][tool]['_'].append(idx)

    with open('metadata/public-server-tools.json', 'w') as handle:
        json.dump(tool_servers, handle)

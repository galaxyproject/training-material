#!/usr/bin/env python
import argparse
import requests
import yaml

from pathlib import Path


class MyDumper(yaml.Dumper):

    def increase_indent(self, flow=False, indentless=False):
        return super(MyDumper, self).increase_indent(flow, False)


def get_metadata_info(tuto, metadata_filepath):
    '''
    Extract the Zenodo DOI of the tutorial from the metadata.yaml file and also 
    the topic title and description

    :param tuto: tutorial name
    :param metadata_filepath: path to the topic metadata.yaml file

    :return: a dictionary with Zenodo record id and DOI, topic title and description
    '''
    meta_info = {'z_record': None, 'z_link': None, 'topic_title': None, 'topic_desc': None}
    with open(metadata_filepath, "r") as m_file:
        metadata = yaml.load(m_file)

        if 'title' not in metadata:
            raise ValueError("No title for the topic in the metadata file")
        meta_info['topic_title'] = metadata['title']

        if 'summary' not in metadata:
            raise ValueError("No summary for the topic in the metadata file")
        meta_info['topic_desc'] = metadata['summary']

        if 'material' not in metadata:
            raise ValueError("No material found in the metadata file")

        for mat in metadata['material']:
            if mat['name'] != tuto:
                continue
            if 'zenodo_link' not in mat or mat['zenodo_link'] == '':
                raise ValueError("Empty Zenodo record found for the tutorial")

            meta_info['z_link'] = mat['zenodo_link']

            if 'doi' in meta_info['z_link']:
                meta_info['z_record'] = meta_info['z_link'].split('.')[-1]
            else:
                meta_info['z_record'] = meta_info['z_link'].split('/')[-1]

    if meta_info['z_record'] is None:
        raise ValueError("No information about the tutorial in the metadata file")

    return meta_info


def create_data_library(meta_info, tuto_dir, overwrite = False):
    '''
    Create the data-library file in the tutorial folder using the information
    on Zenodo (queried though Zenodo API)

    :param meta_info:  a dictionary with Zenodo record id and DOI, topic title and description
    :param tuto_dir: Folder in which the data-library file should be created
    '''
    req = "https://zenodo.org/api/records/%s" % (meta_info['z_record'])
    r = requests.get(req)
    r.raise_for_status()
    req_res = r.json()

    if 'files' not in req_res:
        raise ValueError("No files in the Zenodo record")

    data_lib = {'destination': {'type': 'library',
                                'name': 'GTN - Material',
                                'description': 'Galaxy Training Network Material',
                                'synopsis': 'Galaxy Training Network Material. See https://training.galaxyproject.org'},
                'items': [{'name': '%s' % meta_info['topic_title'],
                          'description': '%s' % meta_info['topic_desc'],
                          'items': []}]}

    for file in req_res['files']:
        file_dict = {'url':'', 'src': 'url', 'ext': '', 'info': meta_info['z_link']}
        if 'type' in file:
            file_dict['ext'] = file['type']
        if 'links' not in file and 'self' not in file['links']:
            raise ValueError("No link for file %s" % file)
        file_dict['url'] = file['links']['self']
        data_lib['items'][0]['items'].append(file_dict)

    data_lib_filepath = tuto_dir / Path("data-library.yaml")
    if data_lib_filepath.is_file():
        print("The data library file already exist and will be overwrite")

    with open(data_lib_filepath, 'w') as stream:
        yaml.dump(data_lib,
                  stream,
                  Dumper=MyDumper,
                  indent=2,
                  default_flow_style=False,
                  default_style='',
                  explicit_start=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create the skeleton of the data-library file for a tutorial with a Zenodo link')
    parser.add_argument('--tutorial', help='Tutorial')
    parser.add_argument('--topic', help='Topic in which the tutorial is')
    #parser.add_argument('--overwrite', help='Topic in which the tutorial is')
    args = parser.parse_args()

    topic_dir = Path("topics") / Path(args.topic)
    if not topic_dir.is_dir():
        raise ValueError("%s is not a topic" % args.topic)
    metadata_filepath = topic_dir / Path("metadata.yaml")
    if not metadata_filepath.is_file():
        raise ValueError("No metadata.yaml file for %s" % args.topic)
    tuto_dir = topic_dir / Path("tutorials") / Path(args.tutorial)
    if not tuto_dir.is_dir():
        raise ValueError("%s is not a tutorial of %s" %(args.tutorial, args.topic))

    meta_info = get_metadata_info(args.tutorial, metadata_filepath)    
    create_data_library(meta_info, tuto_dir)
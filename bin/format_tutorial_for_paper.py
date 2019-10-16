#!/usr/bin/env python

import argparse
import yaml
from pathlib import Path


def check_exists(fp):
    '''Check if file in fp exists
    
    :param fp: Path object to file to test
    '''
    if not fp.exists():
        raise ValueError("%s does not exist" % fp)


def extract_contributors(contributor_fp):
    '''Extract list of contributors
    
    :param fp: Path object to file with contributors
    '''
    contributors = {}
    with open(contributor_fp, 'r') as contributor_f:
        contributors = yaml.load(contributor_f)
    return contributors


def get_author_name(author_id, contributors):
    '''Get author name from list of contributors

    :param author_id: GitHub id of the author
    :param contributors: dictionary with contributor details
    '''
    if author_id in contributors and 'name' in contributors[author_id]:
        return contributors[author_id]['name']
    else:
        return author_id


def get_time(time):
    '''Get time from metadata
    
    :param time: time from metadata
    '''
    if 'h' in time:
        return time.replace('h', ' hours')
    elif 'H' in time:
        return time.replace('H', ' hours')
    elif 'm' in time:
        return time.replace('m', ' minutes')
    elif 'M' in time:
        return time.replace('M', ' minutes')


def create_abstract(metadata):
    '''Create an abstract from the metadata of the tutorial

    :param 
    '''
    a = "<why/aim (1st paragraph of the introduction)>. " \
        "This tutorial is developed for the data analysis platform Galaxy. " \
        "The Galaxy's concept makes high-throughput sequencing data analysis a " \
        "structured, reproducible and transparent process. "
    a += "In this tutorial we focus on %s questions: %s " % (
        len(metadata['questions']), 
        ' '.join(metadata['questions']))
    a += "After finishing this tutorial you will be able to %s." % (
        ', '.join([s[0].lower() + s[1:] for s in metadata['objectives']]))
    a += "Requirements for this tutorial is a minimal Galaxy experience and " \
        "a Galaxy account on the Galaxy Europe server. "
    a += "The estimated duration of this tutorial is around %s." % (
        get_time(metadata['time_estimation']))
        
    return a


def extract_metadata(yaml_content, contributor_fp):
    '''Extract metadata needed by pandoc from tutorial and contributors

    :param yaml_content: YAML content on the top of a tutorial
    :param contributor_fp: Path object to file with contributors
    '''
    metadata = {
        'title': '',
        'author': [],
        'abstract': ''
    }
    contributors = extract_contributors(contributor_fp)
    yaml_metadata = yaml.load(yaml_content)

    metadata['title'] = yaml_metadata['title']
    for c in yaml_metadata['contributors']:
        metadata['author'].append(get_author_name(c, contributors))
    metadata['abstract'] = create_abstract(yaml_metadata)
    
    print(metadata)
    return metadata


def format_tuto_content(content):
    '''Format the tutorial content (boxes)

    :param content: list with tutorial lines
    '''
    return ''.join(content)


def format_tutorial(tuto_fp, formatted_tuto_fp, contributor_fp):
    '''Read tutorial content, format it for pandoc and export it
    
    1. Extract metadata needed by pandoc from tutorial and contributors
    2. Format tutorial content (boxes)
    3. Write in the output file

    :param tuto_fp: Path object to file with original tutorial
    :param formatted_tuto_fp: Path object to file with formatted tutorial for pandoc
    :param contributor_fp: Path object to file with contributors
    '''    
    with open(tuto_fp, "r") as tuto_f:
        full_content = tuto_f.readlines()#yaml.load(tuto_f)

        yaml_content = ''
        for i, l in enumerate(full_content[1:]):
            if l.startswith("---"):
                break
        metadata = extract_metadata(''.join(full_content[1:i]), contributor_fp)
        formatted_tuto_content = format_tuto_content(full_content[i+1:])
        
    with open(formatted_tuto_fp, "w") as tuto_f:
        tuto_f.write('---\n')
        tuto_f.write(yaml.dump(metadata))
        tuto_f.write(formatted_tuto_content)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Format the tutorial file before pandoc')
    parser.add_argument('-t', '--tutorial', help="Path to tutorial directory")
    args = parser.parse_args()

    # get file paths
    tuto_dp = Path(args.tutorial)
    check_exists(tuto_dp)
    original_tuto_fp =  tuto_dp / "tutorial.md"
    check_exists(original_tuto_fp)
    formatted_tuto_fp = tuto_dp / "formatted_tutorial.md"
    contributor_fp = Path("CONTRIBUTORS.yaml")
    check_exists(contributor_fp)

    # 
    format_tutorial(original_tuto_fp, formatted_tuto_fp, contributor_fp)

    
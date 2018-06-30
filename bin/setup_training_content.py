#!/usr/bin/env python
import argparse
import collections
import os
import shutil
from pathlib import Path
import oyaml as yaml


class MyDumper(yaml.Dumper):

    def increase_indent(self, flow=False, indentless=False):
        return super(MyDumper, self).increase_indent(flow, False)


def load_yaml(filepath):
    '''
    Load the content of a YAML file to a dictionary

    :param filepath: path to a YAML file

    :return: dictionary with the content of the YAML file
    '''
    with open(filepath, "r") as m_file:
        content = yaml.load(m_file)
    return content


def save_to_yaml(content, filepath):
    '''
    Save a dictionary to a YAML file

    :param content: dictionary to save
    :param filepath: path an output YAML file
    '''
    with open(filepath, 'w') as stream:
        yaml.dump(content,
                  stream,
                  Dumper=MyDumper,
                  indent=2,
                  default_flow_style=False,
                  default_style='',
                  explicit_start=True)


def change_topic_name(filepath):
    '''
    Change the topic name in the top metadata of a file

    :param filepath: path to index.md, tutorial.md or slides.html
    '''
    with open(filepath, "r") as in_f:
        content = in_f.read()

    content = content.replace("your_topic", args.topic_name)
    content = content.replace("your_tutorial_name", "tutorial1")

    with open(filepath, 'w') as out_f:
        out_f.write(content)


def create_topic(args, topic_dir, template_dir):
    '''
    Create the skeleton of a new topic:

    1. copy templates
    2. update the index.md to match your topic's name
    3. fill the metadata
    4. add a symbolic link to the metadata.yaml from the metadata folder

    :param args: arguments of the script
    :param topic_dir: path to the new topic directory
    :param template_dir: path to the template directory for a new topic
    '''
    # copy templates
    shutil.copytree(template_dir, topic_dir)

    # update the index.md to match your topic's name
    index_path = topic_dir / Path("index.md")
    change_topic_name(index_path)

    # update the metadata file
    metadata_path = topic_dir / Path("metadata.yaml")

    metadata = load_yaml(metadata_path)
    metadata['name'] = args.topic_name
    metadata['title'] = args.topic_title
    metadata['type'] = args.topic_target
    metadata['summary'] = args.topic_summary

    save_to_yaml(metadata, metadata_path)

    # update the metadata in top of tutorial.md and slides.html
    tuto_path = topic_dir / Path("tutorials") / Path("tutorial1")
    hand_on_path = tuto_path / Path("tutorial.md")
    change_topic_name(hand_on_path)
    slides_path = tuto_path / Path("slides.html")
    change_topic_name(slides_path)

    # add a symbolic link to the metadata.yaml
    os.chdir(Path("metadata"))
    os.symlink(Path("..") / metadata_path, Path("{}.yaml".format(args.topic_name)))
    os.chdir(Path(".."))


def update_tuto_file(filepath, keep, args):
    '''
    Update or delete a tutorial (hands-on or slide) file

    :param filepath: path to hands-on or slide file
    :param keep: boolean indicated to keep it or not
    :param args: arguments of the script
    '''
    if keep:
        with open(filepath, "r") as in_f:
            content = in_f.read()

        content = content.replace("your_topic", args.topic_name)
        content = content.replace("your_tutorial_name", args.tutorial_name)

        with open(filepath, 'w') as out_f:
            out_f.write(content)

    elif filepath.is_file():
        filepath.unlink()


def update_tutorial(args, tuto_dir, topic_dir):
    '''
    Update the metadata information of a tutorial

    :param args: arguments of the script
    :param tuto_dir: path to the new tutorial directory
    :param topic_dir: path to the new topic directory
    '''
    # update the metadata file to add the new tutorial
    metadata_path = topic_dir / Path("metadata.yaml")

    metadata = load_yaml(metadata_path)
    found = False
    for mat in metadata["material"]:
        if mat["name"] == args.tutorial_name:
            mat["name"] = args.tutorial_name
            mat["title"] = args.tutorial_title
            mat["hands_on"] = args.tutorial_hands_on
            mat["slides"] = args.tutorial_slides
            found = True

    if not found:
        new_mat = collections.OrderedDict()
        new_mat["title"] = args.tutorial_title
        new_mat["name"] = args.tutorial_name
        new_mat["type"] = 'tutorial'
        new_mat["zenodo_link"] = ''
        new_mat["hands_on"] = args.tutorial_hands_on
        new_mat["slides"] = args.tutorial_slides
        new_mat["workflows"] = False
        new_mat["galaxy_tour"] = False
        new_mat["questions"] = ['', '']
        new_mat["objectives"] = ['', '']
        new_mat["time_estimation"] = '1d/3h/6h'
        new_mat["key_points"] = ['', '']
        new_mat["contributors"] = ['contributor1', 'contributor2']
        metadata["material"].append(new_mat)

    save_to_yaml(metadata, metadata_path)

    # update the metadata in top of tutorial.md or remove it if not needed
    hand_on_path = tuto_dir / Path("tutorial.md")
    update_tuto_file(hand_on_path, args.tutorial_hands_on, args)

    # update the metadata in top of slides.md or remove it if not needed
    slides_path = tuto_dir / Path("slides.html")
    update_tuto_file(slides_path, args.tutorial_slides, args)


def create_tutorial(args, tuto_dir, topic_dir, template_dir):
    '''
    Create the skeleton of a new tutorial

    :param args: arguments of the script
    :param tuto_dir: path to the new tutorial directory
    :param topic_dir: path to the new topic directory
    :param template_dir: path to the template directory for a new tutorial
    '''
    # copy or rename templates
    template_tuto_path = topic_dir / Path("tutorials") / Path("tutorial1")
    if template_tuto_path.is_dir():
        template_tuto_path.rename(tuto_dir)
    else:
        shutil.copytree(template_dir, tuto_dir)

    # fill the metadata of the new tutorial
    update_tutorial(args, tuto_dir, topic_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create the skeleton for a new topic and/or a new tutorial')
    parser.add_argument('--topic_name', help='Name (directory name) of the topic to create or in which the tutorial should be create', required=True)
    parser.add_argument('--topic_title', help='Title of the topic to create', default='Title of the topic', required=False)
    parser.add_argument('--topic_target', help='Target of the topic', default='use', choices=['use', 'admin-dev', 'instructors'], required=False)
    parser.add_argument('--topic_summary', help='Summary of the topic', default='Summary of the topic', required=False)
    parser.add_argument('--tutorial_name', help='Name (directory name) of the new tutorial to create (it will be the directory name)', required=False)
    parser.add_argument('--tutorial_title', help='Title of the new tutorial', default='Title of the tutorial', required=False)
    parser.add_argument('--hands_on', help='Hands-on for the new tutorial', dest='tutorial_hands_on', action='store_true')
    parser.set_defaults(tutorial_hands_on=True)
    parser.add_argument('--slides', help='Slides for the new tutorial', dest='tutorial_slides', action='store_true')
    parser.set_defaults(tutorial_slides=False)
    args = parser.parse_args()

    template_dir = Path("templates")

    topic_dir = Path("topics") / Path(args.topic_name)
    if not topic_dir.is_dir():
        print("The topic {} does not exist. It will be created".format(args.topic_name))
        create_topic(args, topic_dir, template_dir)

    if args.tutorial_name:
        tuto_dir = topic_dir / Path("tutorials") / Path(args.tutorial_name)
        if not tuto_dir.is_dir():
            template_dir = template_dir / Path("tutorials") / Path("tutorial1")
            print("The tutorial {} in topic {} does not exist. It will be created.".format(args.tutorial_name, args.topic_name))
            create_tutorial(args, tuto_dir, topic_dir, template_dir)
        else:
            print("The tutorial {} in topic {} already exists. It will be updated with the other arguments".format(args.tutorial_name, args.topic_name))
            update_tutorial(args, tuto_dir, topic_dir)

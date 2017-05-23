#!/usr/bin/env python

import os
import glob


for topic in os.listdir('./topics'):
    if os.path.isdir( os.path.join('./topics', topic) ):
        topic_path = os.path.join('./topics', topic)
        print(topic_path)

        for dir_name in ['docker', 'slides', 'tutorials']:
            p = os.path.join( topic_path, dir_name)
            assert os.path.exists( p ), '%s not found, but required.' % p
            assert os.path.isdir( p ), '%s not found, but required.' % p

        for file_name in ['metadata.yaml', 'README.md', 'docker/Dockerfile', 'slides/index.html']:
            print(os.path.join( topic_path, file_name))
            p = os.path.join( topic_path, file_name)
            assert os.path.exists( p ), '%s not found, but required.' % p

        for tutorial in os.listdir( os.path.join( topic_path, 'tutorials') ):
            tutorial_path = os.path.join(topic_path, 'tutorials', tutorial)
            print(tutorial)
            if os.path.isdir(tutorial):
                # check for ./topics/<topic>/* files
                for file_name in ["tutorial.md", "metadata.yaml", "tools.yaml", "data-library.yaml"]:
                    p = os.path.join( tutorial_path, file_name)
                    assert os.path.exists( p ), '%s not found, but required.' % p

                # check for ./topics/<topic>/* files
                for dir_name in ["tours", "workflows"]:
                    directory = os.path.join( tutorial_path, dir_name)
                    assert os.path.exists( directory ), '%s does not exists.' % directory
                    assert os.path.isdir( directory ), '%s does not exists.' % directory
                    assert len(os.listdir( directory )) >= 1


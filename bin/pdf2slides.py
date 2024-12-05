#!/usr/bin/env python
# convert pdf slides, into markdown.
#   extract correct images from each page
#   and then the text as well.

import sys
import tempfile
import subprocess
import os
import xml.etree.ElementTree as ET


print("""
---
layout: tutorial_slides

title: ""
questions:
objectives:
key_points:
contributions:
  authorship:
  editing:
---
""")



PDF = os.path.abspath(sys.argv[1])

# pdftohtml -xml -q -stdout CAN_Module1_Lecture.pdf
xml = subprocess.check_output(['pdftohtml', '-xml', '-q', '-stdout', PDF]).decode('utf-8')
tree = ET.fromstring(xml)
for kid in tree:
    sys.stdout.write(f"\n\n--- # {kid.attrib['number']}\n\n")
    lasttop = 0
    for e in kid:
        if e.tag == 'fontspec': continue

        if e.tag == "image":
            print(f"![](images/{os.path.basename(e.attrib['src'])})")
        elif e.tag == 'text':
            if e.text is not None:
                if abs(int(e.attrib['top']) - lasttop) < 5:
                    sys.stdout.write(" " + e.text)
                else:
                    sys.stdout.write("\n" + e.text)
            else:
                if len(e) == 1:
                    if abs(int(e.attrib['top']) - lasttop) < 5:
                        sys.stdout.write(" " + ET.tostring(e[0]).decode('utf-8'))
                    else:
                        sys.stdout.write("\n" + ET.tostring(e[0]).decode('utf-8'))
                else:
                    raise Exception("Unsupported")
        else:
            raise Exception(f"Unsupported tag {e.tag}")
        lasttop = int(e.attrib['top'])

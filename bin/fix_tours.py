import os
import sys
import yaml
from bioblend import toolshed
from urllib.parse import unquote as decode

'''
A program to update outdated Galaxy Tours to conform to the latest version 
(21.01 as of this writing) of the Galaxy UI.

This program generates a Jinja2 template and externalizes tool IDs to a 
separate vars.yaml file so tool versions can be updated automatically (see the
update-tour.py program) without editing the tour file itself.

A list of the tools used in the tour are  written to a tools.yaml file so 
admins can easily see and install the tools needed by a particular tour. The
tools.yaml file is also updated by the update-tour.py program.

Use the render-tour.py program to generate the tour YAML file that will be 
loaded by Galaxy.

NOTE: Since YAML is a PITA to roundtrip cleanly and easily we don't use a YAML
parser to read the original tour.yaml file, but simple read it line by line and
do string replacements.  The problem is Block Literals, i.e. multi-line strings
specified with >- or |-. Once those strings are stored in dictionaries the style
of input is lost and the YAML processor will pick what it thinks is best when
it outputs the string.
'''

DEFAULT_TOOLSHED = 'https://toolshed.g2.bx.psu.edu'
shed = toolshed.ToolShedInstance(DEFAULT_TOOLSHED)

# The strings that we do simple search/replace on.
# TODO These values should be parameterized.
replacements = [
    (' .fa.fa-upload', ''),
    ('#history-options-button', '#history-new-button'),
    ('#tool-search-query', '.search-query')
]


def make_tool_entry(tool_id):
    '''Create a dictionary to generated the YAML for the tools.yaml file.'''
    id = decode(tool_id)
    parts = id.split('/')
    if len(parts) == 1:
        return {
            'name': id,
            'owner': None,
            'tool_shed_url': DEFAULT_TOOLSHED,
            'tool_panel_section_label': 'Tours'
        }
    return {
        'name': parts[3],
        'owner': parts[2],
        'tool_shed_url': 'https://' + parts[0],
        'tool_panel_section_label': 'Tours'
    }


def parse_name_and_tag(file_path):
    '''Extract a name and tag from the file path.'''
    parts = file_path.split('/')
    return (parts[1], parts[3])


def do_replacement(s):
    for before, after in replacements:
        s = s.replace(before, after)
    return s


def get_indentation(line):
    index = len(line) - len(line.lstrip())
    return line[:index]


def parse_tour(tour_file, directory):
    print(f'Parsing {tour_file}')
    vars = {}
    lines = []
    tools = []
    seen = []
    tag, name = parse_name_and_tag(tour_file)
    with open(tour_file) as f:
        for line in f.readlines():
            line = do_replacement(line)
            if line.startswith("steps:"):
                lines.append("tags:\n")
                lines.append(f'  - "{ tag }"\n')
                lines.append('  - auto\n')
                lines.append(line)
            elif line.find("/tool_runner?tool_id=") > 0:
                start = line.index("tool_id=") + 8
                end = line.index('"', start)
                id = line[start:end]
                indent = get_indentation(line)
                # id = decode(click[start:end])
                tool_info = make_tool_entry(id)
                var_name = tool_info['name']
                if '+' in var_name:
                    var_name = var_name.replace('+','_')
                new_line = indent + 'a[href$="/tool_runner?tool_id={{ ' + var_name + ' }}"]'
                if line.endswith('.tool-old-link'):
                    new_line += ' .tool-old-link'
                lines.append(new_line + '\n')
                if tool_info not in tools:
                    tools.append(tool_info)
                vars[var_name] = id
            else:
                lines.append(line)
    lines.append('\n')
    output_path = os.path.join(directory, name)
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    with open(f"{output_path}/tour.yaml.j2", "w") as f:
        f.writelines(lines)
    with open(f"{output_path}/vars.yaml", "w") as f:
        yaml.dump(vars, f)
    with open(f"{output_path}/tools.yaml", "w") as f:
        yaml.dump({ 'tools':tools }, f)


def fix_tours(directory):
    for root, dirs, files in os.walk('topics'):
        for file in files:
            if file == 'tour.yaml':
                parse_tour(os.path.join(root, file), directory)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"USAGE: {sys.argv[0]} /output/directory")
        sys.exit(1)
    directory = sys.argv[1]
    if not os.path.exists(directory):
        print(f'Could not find {directory}')
        sys.exit(1)
    if not os.path.isdir(directory):
        print('The specified path must be a directory')
        sys.exit(1)
    fix_tours(directory)

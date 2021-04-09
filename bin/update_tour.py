import os
import sys
import yaml
from bioblend import toolshed
from urllib.parse import quote_plus as urlencode

'''
Update the tool versions used in a given tour.  Both the tools.yaml and
vars.yaml files are updated with the latest version available on the ToolShed.
'''

DEFAULT_TOOLSHED = 'https://toolshed.g2.bx.psu.edu'

# Common keys into the tools dict. Defined solely so our IDE can do completions
# and I don't consistently misspell revisisions or have to remember if it is
# toolshed_url or tool_shed_url
NAME = 'name'
OWNER = 'owner'
TOOLS = 'tools'
SHED = 'tool_shed_url'
REVISIONS = 'revisions'

tool_sheds = {DEFAULT_TOOLSHED: toolshed.ToolShedInstance(DEFAULT_TOOLSHED)}

def get_guid(info):
    '''
    Get the GUID for a given tool.

    The tool's GUID is the tool id used by selectors in the tour.
    '''
    for m in info:
        if 'model_class' in m and m['model_class'] == 'RepositoryMetadata':
            return m['valid_tools'][0]['guid']


def validate(tool):
    """Ensure the tool has the fields we need so we don't need to check later."""
    if SHED not in tool:
        tool[SHED] = DEFAULT_TOOLSHED
    if REVISIONS not in tool:
        tool[REVISIONS] = []


def get_tool_shed(tool):
    '''
    Get the ToolShed URL for the toolshed the tool was installed from.

    #TODO We shoud probably raise an exception if a tool comes from anywhere
          but the main tool shed.
    '''
    url = tool[SHED]
    if url in tool_sheds:
        ts = tool_sheds[url]
    else:
        ts = toolshed.ToolShedInstance(url)
        tool_sheds[url] = ts
    return ts


def update_tour(tour_dir):
    '''
    Update the tools.yaml and vars.yaml with the latest tool versions used
    in the tour.
    '''
    print(f'Updating {os.path.basename(tour_dir)}')
    tools_file = os.path.join(tour_dir, 'tools.yaml')
    if not os.path.exists(tools_file):
        print(f'Could not find the tools.yaml file.')
        return
    vars_file = os.path.join(tour_dir, 'vars.yaml')
    if not os.path.exists(vars_file):
        print(f'Could not find the vars.yaml file')
        return
    with open(vars_file, 'r') as f:
        variables = yaml.safe_load(f)
    with open(tools_file, "r") as f:
        data = yaml.safe_load(f)
    tool_list = data[TOOLS]
    print(f'Updating {len(tool_list)} tools.')
    for tool in tool_list:
        name = tool[NAME]
        owner = tool[OWNER]
        print(f"Getting latest revision for {name}")
        validate(tool)
        ts = get_tool_shed(tool)
        revs = ts.repositories.get_ordered_installable_revisions(name, owner)
        if revs and len(revs) > 0:
            tool[REVISIONS] = [ revs[-1] ]
            info = ts.repositories.get_repository_revision_install_info(name, owner, revs[-1])
            guid = get_guid(info)
            variables[name] = urlencode(guid)

    data = {"tools": tool_list}
    with open(tools_file, "w") as f:
        yaml.dump(data, f)
    with open(vars_file, "w") as f:
        yaml.dump(variables, f)
    print(f"Tour {tour_dir} updated")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f'USAGE: python {sys.argv[0]} /path/to/tour_dir/')
        sys.exit(1)

    path = sys.argv[1]
    if not os.path.exists(path):
        print(f'The path "{path}" does not appear to contain a Galaxy Tour.')
        sys.exit(1)
    update_tour(path)

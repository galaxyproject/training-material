import os
import yaml
import sys

HEAD = """
---
api_key: admin
galaxy_instance: http://localhost:8080
install_resolver_dependencies: true
install_tool_dependencies: false
"""
tutorial = str(sys.argv[1])
workflow = str(sys.argv[2])
workflowdir = tutorial + "/workflows/"
filename = workflowdir + "wftools.yaml"

os.system("workflow-to-tools -w {} -o {} -l {}".format(workflowdir + workflow, filename, os.path.basename(os.path.normpath(tutorial))))

yamltoolf = yaml.load(open(filename))
with open(filename, 'w') as toolfile:
    toolfile.write(HEAD)
    toolfile.write(yaml.dump({'tools':yamltoolf['tools']}, default_flow_style=False, default_style=''))
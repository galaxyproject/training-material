import os
import sys
import yaml
from jinja2 import Template


def get_tour_id(template_file):
    dirname = os.path.dirname(template_file)
    return os.path.basename(dirname)


def render(template_file):
    dir = os.path.dirname(template_file)
    vars_file = os.path.join(dir, 'vars.yaml')
    with open(template_file) as f:
        template = Template(f.read())
    with open(vars_file) as f:
        vars = yaml.safe_load(f)
    # id = get_tour_id(template_file)
    output = template.render(vars)
    outfile = os.path.join(dir, f"tour.yaml")
    with open(outfile, 'w') as f:
        f.write(output)
    print(f'Rendered {outfile}')


if __name__ == '__main__':
    if len(sys.argv) == 1 or len(sys.argv) > 3:
        print(f'USAGE: python {sys.argv[0]} /path/to/tour.yaml.j2')
        print('\nThe rendered tour will be written to the same directory as the template file.\n')
        sys.exit(1)

    path = sys.argv[1]
    if not path.endswith('tour.yaml.j2'):
        if os.path.isdir(path):
            path = os.path.join(path, 'tour.yaml.j2')
        else:
            print(f'The path "{path}" does not appear to contain Galaxy Tour.')
            sys.exit(1)
    if not os.path.exists(path):
        print(f'The path "{path}" does not appear to contain a Galaxy Tour.')
        sys.exit(1)
    render(path)
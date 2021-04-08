import os
import sys

'''
At least one tour (iwtomics) uses incorrect indentation that causes yamllint
to complain.
'''

def fix(path):
    with open(path) as f:
        for line in f.readlines():
            if len(line) > 4 and line.startswith('  '):
                print(line[2:], end='')
            else:
                print(line, end='')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('No filename provided.')
        sys.exit(1)
    fix(sys.argv[1])
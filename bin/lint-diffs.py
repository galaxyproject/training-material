#!/usr/bin/env python3
import sys
import argparse
parser = argparse.ArgumentParser(description='Check for missing whitespace on specific lines.')
parser.add_argument('tutorial', type=argparse.FileType('r'))
parser.add_argument('--fix', action='store_true')
args = parser.parse_args()


def lint(fh):
    in_diff = False
    min_length = sys.maxsize
    failures = 0
    fixed = []
    for lineno, line in enumerate(fh.readlines()):
        line = line.rstrip('\n')
        if '```diff' in line:
            in_diff = True
            min_length = line.index('`') + 1
        elif in_diff and '```' in line:
            in_diff = False
            min_length = sys.maxsize
        elif in_diff and len(line) < min_length:
            print(f'{fh.name}: {lineno}: diff line too short ({len(line)} < {min_length})!: "{line}"')
            line += " " * (min_length - len(line))
            failures += 1
        fixed.append(line)
    return failures, fixed


if __name__ == '__main__':
    failures, fixed = lint(args.tutorial)
    if args.fix:
        with open(args.tutorial.name, 'w') as handle:
            handle.write("\n".join(fixed))
    if failures:
        print(f'ERROR: Linting diffs failed with {failures} failures')
    print('All diffs OK!')

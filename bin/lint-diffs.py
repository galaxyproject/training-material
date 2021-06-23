#!/usr/bin/env python3
import sys
import fileinput


def lint(fh):
    in_diff = False
    min_length = sys.maxsize
    failures = 0
    for line in fh:
        line = line.rstrip('\n')
        if '```diff' in line:
            in_diff = True
            min_length = line.index('`') + 1
        elif in_diff and '```' in line:
            in_diff = False
            min_length = sys.maxsize
        elif in_diff and len(line) < min_length:
            print(f'{fh.filename()}: {fileinput.filelineno()}: diff line too short ({len(line)} < {min_length})!: "{line}"')
            failures += 1
    return failures
        
            
def main(): 
    failures = None
    with fileinput.input() as fh:
        failures = lint(fh)
    if failures:
        print(f'ERROR: Linting diffs failed with {failures} failures')
        sys.exit(1)
    print('All diffs OK!')
    sys.exit(0)


main()   

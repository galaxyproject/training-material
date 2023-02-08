#!/bin/bash

# GTN AutoFixer script. This fixes any automatically correctable errors we can
# detect.

# faqs/Must be called with 'snippet', not 'include'
sed -i "s/{% include faqs/{% snippet faqs/g" $(grep '{% include faqs' -R topics -l)

# gat-diffs/These can be auto-corrected it seems
find ./topics/admin/ -name '*.md' -type f -print0 | xargs -n 1 -0 python3 bin/lint-diffs.py

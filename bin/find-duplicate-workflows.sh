#!/bin/bash
md5sum $(find topics -name '*.ga') | grep $(md5sum $(find topics -name '*.ga') | sort | cut -f 1 -d' ' | sort | uniq -c | sort -n | grep -v '^\s*1 ' | awk '{print $2}')

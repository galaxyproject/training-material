#!/bin/bash
grep '(missing reference)' -R _site && echo "Missing references detected" && exit 1

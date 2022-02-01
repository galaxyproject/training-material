#!/usr/bin/env python
import sys
import yaml
import json

[sys.stdout.write(json.dumps(doc, indent=2)) for doc in yaml.safe_load_all(sys.stdin)]

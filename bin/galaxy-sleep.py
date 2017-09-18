#!/usr/bin/env python
import sys
import requests
import time
count = 0

while True:
    try:
        result = requests.get(sys.argv[1] + '/api/version').json()
        sys.stdout.write("Galaxy Version: %s\n" % result['version_major'])
        break
    except Exception as e:
        sys.stdout.write("[%02d] Galaxy not up yet... %s\n" % (count, str(e)[0:100]))
        sys.stdout.flush()

    count += 1
    time.sleep(1)

sys.exit(0)

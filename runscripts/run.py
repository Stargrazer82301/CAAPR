import os
import sys
import json
import CAAPR

if len(sys.argv) <= 1:  # No config file given as argument: use default
    config_filename = './runconfig.json'
else:
    config_filename = sys.argv[1]
if not os.path.exists(config_filename):
    raise ValueError("Can not find config file at " + config_filename)

with open(config_filename, 'r') as config_file:
    d_config = json.load(config_file)

print('Starting CAAPR with following config: ' + str(d_config))
CAAPR.CAAPR_Main.Run(**d_config)

# Jubilate
print('All done!')
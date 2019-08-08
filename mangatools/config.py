#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

ESP = 1e-8 # the error

release = 'MPL-8'
print("Data relese: {}, configure: mangatools/config.py".format(release))

# Data release MPL-7
if release == 'MPL-7':
    DRP_VERSION = 'v2_4_3' # the default drp version
    DAP_VERSION = '2.2.1'
elif release == 'MPL-8':
# MPL-8
    DRP_VERSION = 'v2_5_3'
    DAP_VERSION = '2.3.0'
PRODOCTS = 'HYB10'

## Data directory, you can direct change specific path after import this module
SAS = os.getenv('SAS_BASE_DIR', default=os.path.expanduser('~')+'/SAS')

print("Global SAS directory is {0}".format(SAS))

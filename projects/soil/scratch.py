#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:21:14 2020

@author: travis
"""

from functions import get_gssurgo, get_gnatsgo, map_variable
from functions import mukey_rescue

# Sample arguments for map_variable
MUKEY = "~/data/weto/soil/mukey_de.tif"
GDB = "~/data/weto/soil/gNATSGO_DE.gdb"
DST = "~/data/weto/soil/brockdepmin_de.tif"
VARIABLE = "brockdepmin"


if not os.path.exists(os.path.expanduser(MUKEY)):
    get_gssurgo()
    get_gnatsgo()

mukey_rescue(GDB, "0", MUKEY)
map_variable(gdb_path=GDB,
             mukey_path=MUKEY,
             variable=VARIABLE,
             dst=DST)



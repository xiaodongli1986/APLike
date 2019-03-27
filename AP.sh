#!/bin/bash
APPATH=/home/xiaodongli/software/APLike

# libarary
export LIBRARY_PATH=$APPATH/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=$APPATH/lib:$LD_LIBRARY_PATH

# moduls path
export APmods=$APPATH/mods

# link libarary, include modules
export APlm=-lAP\ -I$APmods


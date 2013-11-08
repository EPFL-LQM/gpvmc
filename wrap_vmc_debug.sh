#!/bin/bash
# For allowing core dumps on the remote hosts to be
# produced if a program aborts or segfaults
ulimit -c unlimited
./vmc "$@"

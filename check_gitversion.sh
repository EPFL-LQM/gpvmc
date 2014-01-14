#!/bin/bash

gv=`git rev-parse --short HEAD`
cv=""
if [ -e gitversion.h ]
then
    cv=`grep -E gitversion=\"[[:alnum:]]+\" gitversion.h | cut -d \" -f 2`
fi
if [ "$gv" != "$cv" ]
then
    sed "s/\".*\"/\"`git rev-parse --short HEAD`\"/" <gitversion.h.tpl >gitversion.h
fi

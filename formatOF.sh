#!/bin/sh
find . -iname *.C -o -iname *.H | xargs clang-format -i --style=file:"$1"
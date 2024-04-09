#!/bin/sh
find applications -iname *.C -o -iname *.H | xargs clang-format -i --style=file:"$1"
find etc -iname *.C -o -iname *.H | xargs clang-format -i --style=file:"$1"
find modules -iname *.C -o -iname *.H | xargs clang-format -i --style=file:"$1"
find src -iname *.C -o -iname *.H | xargs clang-format -i --style=file:"$1"
find tutorials -iname *.C -o -iname *.H | xargs clang-format -i --style=file:"$1"
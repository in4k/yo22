#!/bin/sh

cc -Werror -Wall -O0 -g -lX11 -lGL -DTOOL=1 -DDEBUG -lm yo22.c -o yo22_tool && $1 ./yo22_tool

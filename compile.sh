#!/bin/sh
set -eu
cc -O2 -std=c11 -Wall -Wextra spookyLights3D.c scene.c utils.c svdDynamic.c -lm -o spookyLights3D

#!/usr/bin/env bash


mkdir Refs
cd Refs

curl -L -o out.zip "https://www.dropbox.com/sh/y5ly6bswf4qqrg9/AACWrqJ74AN8RRPLVD9xe38Fa?dl=1" --http1.1

unzip out.zip


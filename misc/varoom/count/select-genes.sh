#!/bin/bash

# Take an input bed file, and pick a random gene from each chromosome

sort -R $1 | awk '/[-,_]/ { next } {a[$1] = $4} END {for (k in a) {print a[k]}}' | sort

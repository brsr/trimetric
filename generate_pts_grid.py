#!/usr/bin/env python3
# -*- coding: utf-8 -*-
stepsize = 0.1
for x in range(int(-180/stepsize), int(180/stepsize) + 1):
    for y in range(int(-90/stepsize), int(90/stepsize) + 1):
        print("{0:.1f} {1:.1f}".format(x*stepsize, y*stepsize))
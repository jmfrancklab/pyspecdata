#!/bin/bash
while [ 1 ]; do
    sleep 1s
    make -k notebook_wc.pdf # | tail
done

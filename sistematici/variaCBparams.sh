#!/bin/sh

for ((i=0; i<10; i++))
do
    python variaCBparams.py -f --makeTree $i
done


for ((i=0; i<10; i++))
do
    python variaCBparams.py -f --makeTree $i -b1
done


for ((i=0; i<10; i++))
do
    python variaCBparams.py -f --makeTree $i -b2
done


for ((i=0; i<10; i++))
do
    python variaCBparams.py -f --makeTree $i -b3
done


for ((i=0; i<10; i++))
do
    python variaCBparams.py -f --makeTree $i -b4
done



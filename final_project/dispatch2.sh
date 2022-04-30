#!/bin/sh

host="simpson"
num=21

while [ $num -lt 29 ]
do
  scp -o StrictHostKeyChecking=no $host$num:blackhole_render*.xwd .
  num=`expr $num + 1`
done

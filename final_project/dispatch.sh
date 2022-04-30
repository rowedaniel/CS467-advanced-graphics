#!/bin/sh

host="simpson"
num=21

while [ $num -lt 30 ]
do
  scp -o StrictHostKeyChecking=no "./a.out" $host$num:~/OLDF21/deleteme.out
  num=`expr $num + 1`
done

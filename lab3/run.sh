#!/bin/bash
N1=1
N2=60000
delta=600001
for ((N=$N1; N<=$N2; N=(($N+$delta)))) do
    echo $(./lab3.exe $N 1) >> out/lab3-1.txt
    echo $(./lab3.exe $N 2) >> out/lab3-4.txt
    echo $(./lab3.exe $N 4) >> out/lab3-8.txt
    echo $(./lab3.exe $N 6) >> out/lab3-12.txt
done


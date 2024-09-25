N1=1500
N2=61000
delta=5950

for ((N=$N1; N<=$N2; N=(($N+$delta)))) do
    echo $(./lab3-seq $N) >> out/lab3-seq.txt
    echo $(./lab3-par-1 $N) >> out/lab3-1.txt
    echo $(./lab3-par-10 $N) >> out/lab3-10.txt
    echo $(./lab3-par-12 $N) >> out/lab3-12.txt
    echo $(./lab3-par-20 $N) >> out/lab3-20.txt
done
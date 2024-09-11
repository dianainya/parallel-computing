N1=1500
N2=61000
delta=5950

for ((N=$N1; N<=$N2; N=(($N+$delta)))) do
    echo $(./lab1-seq $N) >> out/lab1-seq.txt
    echo $(./lab1-par-1 $N) >> out/lab1-1.txt
    echo $(./lab1-par-10 $N) >> out/lab1-10.txt
    echo $(./lab1-par-12 $N) >> out/lab1-12.txt
    echo $(./lab1-par-20 $N) >> out/lab1-20.txt
done
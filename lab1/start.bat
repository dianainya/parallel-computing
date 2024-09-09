@echo off
setlocal enabledelayedexpansion

set N1=1000
set N2=32200
set delta=3120

for /l %%N in (%N1%,1,%N2%) do (
    set /a mod=%%N %% %delta%
    if !mod! EQU 0 (
        lab1-seq.exe %%N >> out\lab1-seq.txt
        lab1-par-1.exe %%N >> out\lab1-1.txt
        lab1-par-10.exe %%N >> out\lab1-10.txt
        lab1-par-12.exe %%N >> out\lab1-12.txt
        lab1-par-20.exe %%N >> out\lab1-20.txt
    )
)

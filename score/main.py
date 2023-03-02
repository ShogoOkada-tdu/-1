import sys
import fileinput
from pathlib import Path
import pandas as pd
import numpy as np
import re
import math


#ファイルの読み込み

args = sys.argv
print ("Motif : " + args[1])

if Path(sys.argv[1]).exists():
        with open(args[1], 'r') as f:
            motif = f.read().split("\n")
            #print(motif)

else:  
    print("not exists")

if Path(sys.argv[2]).exists():
    with open(args[2], 'r') as f:
        DNA = f.read().split("\n")
else:  
    print("")


#頻度表の作成

print('頻度表')

A = []
C = []
G = []
T = []

for i in range(len(motif)):
    for j in motif[i]:
        if j == 'A':
                A.append(1)
                C.append(0)
                G.append(0)
                T.append(0)
        elif j == 'C':
                A.append(0)
                C.append(1)
                G.append(0)
                T.append(0)
        elif j == 'G':
                A.append(0)
                C.append(0)
                G.append(1)
                T.append(0)
        elif j == 'T':
                A.append(0)
                C.append(0)
                G.append(0)
                T.append(1)
            

print('Aの個数')
length = len(A)
n = 0
s = len(motif[0])
a = [0] * len(motif[0])
for i in a:
    a = a + np.array(A[n:n + s :1])
    n += s
    if n >= length:
        print(a)
        break

print('Cの個数')
length = len(C)
n = 0
s = len(motif[0])
c = [0] * len(motif[0])
for i in C:
    c = c + np.array(C[n:n + s :1])
    n += s
    if n >= length:
        print(c)
        break
 

print('Gの個数')
length = len(G)
n = 0
s = len(motif[0])
g = [0] * len(motif[0])
for i in G:
    g = g + np.array(G[n:n + s :1])
    n += s
    if n >= length:
        print(g)
        break

print('Tの個数')
length = len(T)
n = 0
s = len(motif[0])
t = [0] * len(motif[0])
for i in T:
    t = t + np.array(T[n:n + s :1])
    n += s
    if n >= length:
        print(t)
        print('')
        break


#scoreの計算

print('scoreの計算')

Ag = 7519429
Cg = 4637676
Gg = 4637676
Tg = 7519429
total = 24314210

As = []
Cs = []
Gs = []
Ts = []

for i in range(len(motif[0])):
    Ak = np.array(a[i]) + 1
    Ab = np.array(a[i]) + np.array(c[i]) + np.array(g[i]) + np.array(t[i]) + 4
    Akb = Ak / Ab
    ss = Ag / total
    log = math.log(Akb / ss)
    As.append(round(log, 2))

print(As)

for i in range(len(motif[0])):
    Ck = np.array(c[i]) + 1
    Cb = np.array(a[i]) + np.array(c[i]) + np.array(g[i]) + np.array(t[i]) + 4
    Ckb = Ck / Cb
    ss = Cg / total
    log = math.log(Ckb / ss)
    Cs.append(round(log, 2))

print(Cs)

for i in range(len(motif[0])):
    Gk = np.array(g[i]) + 1
    Gb = np.array(a[i]) + np.array(c[i]) + np.array(g[i]) + np.array(t[i]) + 4
    Gkb = Gk / Gb
    ss = Gg / total
    log = math.log(Gkb / ss)
    Gs.append(round(log, 2))

print(Gs)

for i in range(len(motif[0])):
    Tk = np.array(t[i]) + 1
    Tb = np.array(a[i]) + np.array(c[i]) + np.array(g[i]) + np.array(t[i]) + 4
    Tkb = Tk / Tb
    ss = Tg / total
    log = math.log(Tkb / ss)
    Ts.append(round(log, 2))

print(Ts)
print()


#結合部位の探索
print('結合部位の探索')
hit = 0
for i in range(len(motif)):
    for j in range(len(motif[0])):
        if motif[i][j] == 'A':
            hit = hit + As[j]
        elif motif[i][j]  == 'C':
            hit = hit + Cs[j]
        elif motif[i][j]  == 'G':
            hit = hit + Gs[j]
        elif motif[i][j]  == 'T':
            hit = hit + Ts[j]
    print(motif[i])
    print(round(hit, 2))
    hit = 0






    
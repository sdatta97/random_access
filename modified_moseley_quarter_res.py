import math
import numpy as np

V0 = [0.0, 0.0]
V1 = [[0.0] * 385 for _ in range(2)]
V2 = [[[0.0] * 589 for _ in range(589)] for _ in range(2)]
V11 = [[[0.0] * 385 for _ in range(385)] for _ in range(2)]
V12 = [[[0.0] * 385 for _ in range(385)] for _ in range(2)]
alph = [0.0] * 385
beta = [[0.0] * 589 for _ in range(589)]
NUM_ITERATIONS = 50
vmax = 0.0
tmax1 = 0.99
tex = 0.0
t = 0.0
tmax = 0.0
vn = 0.0
vn1 = 0.0
vn2 = 0.0
tdif = 0.0
tdif1 = 0.0
aux = 0.0
xmin = 0.0
xmax = 0.0
x = 0.0
dif = 0.0
dif1 = 0.0
b1max = 0.0
bxdif = 0.0
bxdif1 = 0.0
exy = 0.0
bmax = 0.0
bmax1 = 0.0
xy = 0.0
ex = 0.0
exb = 0.0
exyb = 0.0
axdif = 0.0
axdif1 = 0.0
a = 0.0
ea = 0.0
ea2 = 0.0
eax = 0.0
eax2 = 0.0
amax1 = 0.0
amax = 0.0
adif = 0.0
adif1 = 0.0
b = 0.0
y = 0.0
ey = 0.0
eb = 0.0
eyb = 0.0
bdif = 0.0
bdif1 = 0.0
i = 0
j = 0
k = 0
l = 0
m = 0
n = 0
it = 0
ib = 0
iby = 0
lb = 0
ia = 0
iax = 0
la = 0
k1 = 0
kb = 0
ibxy = 0
ibx = 0
kl = 0
m1 = 1
m2 = m1 + 1
del_ = 0.0025
delf = 0.00025
del1 = del_ * m1
delf1 = delf * m1
flag = 0
type_ = 0
V1[0][0] = 0.0
V0[1] = 0.0
for j in range(385):
    V1[1][j] = 0.0
for j in range(589):
    for k in range(589):
        V2[1][j][k] = 0.0
for j in range(385):
    for k in range(385):
        V12[1][j][k] = 0.0
        V11[1][j][k] = 0.0
for n in range(1, NUM_ITERATIONS + 1):
    i = 1 - (n - (n // 2) * 2)
    j = 0 + n - (n // 2) * 2
    print("i: %f, j: %f" %(i,j))
    vmax = 0.0
    tmax1 = 0.99
    for l in range(1, 189):
        t = tmax1 + l * del_
        tex = math.exp(-t)
        it = 400 + l
        vn = t * tex + (1 + t) * tex * V0[j] + V2[j][it - 1][0] * (1 - (1 + t) * tex)
        if vn > vmax:
            vmax = vn
            tmax = t
    print("Course S0: %f, %f" % (vmax, tmax))
    for l in range(1, 20):
        t = tmax - del_ + l * delf
        if t < 0:
            continue
        tex = math.exp(-t)
        aux = t / del_
        it = int(np.floor(aux))
        dif = aux - it
        dif1 = 1 - dif
        it = it + 1
        vn = t * tex + (1 + t) * tex * V0[j] + (V2[j][it - 1][0] * dif1 + V2[j][it][0] * dif) * (1 - (1 + t) * tex)
        if vn > vmax:
            tmax = t
            vmax = vn
    V0[i] = vmax
    t = tmax
    print("Fine S0: %f, %f" % (vmax, tmax))
    x = V0[i] - V0[j]
    xmin = x
    xmax = x
    print("Value change S0: %f" % x)
    for l in range(m2, 385, m1):
        y = (l - 1) * del_
        ey = math.exp(-y)
        vn = 1 + V1[j][l - 1]
        b = 0.0
        bmax = 0.0
        vmax = vn
        for b in np.arange(delf1, y, delf1):
            eb = math.exp(-b)
            eyb = ey / eb
            aux = b / del1
            ib = int(np.floor(aux))
            dif = aux - ib
            dif1 = 1 - dif
            ib = ib * m1 + 1
            vn = (eb - ey) * (1 + V1[j][ib - 1] * dif + V1[j][ib - 1 + m1] * dif1) + (1 - eb) * V2[j][0][ib - 1 + m1] * dif
            if ib == 1:
                vn = vn + (eyb - ey) * 0.5 * (1 + V1[j][0] + V2[j][0][0]) * dif1
            if ib > 1:
                vn = (vn + (1 - eb) * V2[j][0][ib - 1] * bdif1) / (1 - ey)
            if vn > vmax:
                vmax = vn
                bmax = b
        beta[0][l - 1] = bmax
        V2[i][0][l - 1] = vmax
        x = V2[i][0][l - 1] - V2[j][0][l - 1]
        if x < xmin:
            xmin = x
        if x > xmax:
            xmax = x
        if xmin < 0:
            print("Neg xmin: %f" % xmin)
    vmax = 0.0
    for l in range(1, 386, m1):
        a = (l - 1) * del_
        if a < 0:
            continue
        ea = math.exp(-a)
        vn = ea * (1 + V0[i]) + a * ea * (1 + V1[j][0]) + (1 - ea - a * ea) * (V12[j][0][l - 1])
        if vn > vmax:
            amax = a
            vmax = vn
    for l in range(1, 20):
        a = amax - del1 + l * delf1
        if a < 0:
            continue
        ea = math.exp(-a)
        aux = a / del1
        ia = int(np.floor(aux))
        dif = aux - ia
        dif1 = 1 - dif
        ia = ia * m1 + 1
        if (ia < 385):
            vn = ea * (1 + V0[i]) + a * ea * (1 + V1[j][0]) + (1 - ea - a * ea) * (V12[j][0][ia + m1 - 1] * adif + V12[j][0][ia - 1] * adif1)
        else:
            vn = ea * (1 + V0[i]) + a * ea * (1 + V1[j][0]) + (1 - ea - a * ea) * V12[j][0][ia]
        if vn > vmax:
            amax = a
            vmax = vn
    alph[0] = amax
    V1[i][0] = vmax
    x = V1[i][0] - V1[j][0]
    if x < xmin:
        xmin = x
    if x > xmax:
        xmax = x
    if xmin < 0:
        print("Neg xmin: %f" % xmin)
    for la in range(m2, 386, m1):
        a = (la - 1) * del_
        ea = math.exp(-a)
        ea2 = math.exp(-0.5 * a)
        V12[i][0][la - 1] = (ea2 * (1 - (1 + a / 2) * ea2) * V12[j][0][(la - 1) // 2] + 0.5 * a * ea2 * (1 - ea2) * V11[j][0][(la - 1) // 2] + (1 - (1 + a / 2) * ea2) * V12[j][0][(la - 1) // 2]) / (1 - (1 + a) * ea)
        V11[i][0][la - 1] = (a * ea * V1[j][0] + (1 - (1 + a) * ea) * V12[j][0][la - 1]) / (1 - ea)
    for l in range(m2, 386, m1):
        x = (l - 1) * del_
        ex = math.exp(-x)
        V12[i][l - 1][0] = 0.5 * (V12[j][l - 1][0] + V11[j][l - 1][0])
        V11[i][l - 1][0] = 0.5 * (V1[j][l - 1] + V12[j][l - 1][0])
    V12[i][0][0] = 0.5 * (V12[j][0][0] + V11[j][0][0])
    V11[i][0][0] = 0.5 * (V1[j][0] + V12[j][0][0])
    V2[i][0][0] = (V2[j][0][0] + V1[j][0] + 1) / 2
    x = V2[i][0][0] - V2[j][0][0]
    if x < xmin:
        xmin = x
    if x > xmax:
        xmax = x
    if xmin < 0:
        print("Neg xmin: %f" % xmin)
    for k in range(m2, 386, m1):
        x = (k - 1) * del_
        ex = math.exp(-x)
        for l in range(1, 386 - k, m1):
            y = (l - 1) * del_
            ey = math.exp(-y)
            exy = ex * ey
            b = 0.0
            vn = 0.0
            bmax = 0.0
            b1max = 0.0
            vmax = vn
            flag = 0
            type_ = 0
            for b in np.arange(delf1, x + y, delf1):
                if b >= x:
                    flag = 1
                    break
                eb = math.exp(-b)
                exb = ex / eb
                exyb = ex * ey / eb
                aux = b / del1
                ib = int(np.floor(aux))
                dif = aux - ib
                dif1 = 1 - dif
                ib = ib * m1 + 1
                ibxy = k + l - ib - m1
                ibx = k - ib - m1 + 1
                vn = (eb - ex - (x - b) * exy) * (V2[j][ibx - 1][l - 1] * dif + V2[j][ibx + m1 - 1][l - 1] * dif1) + b * (eb - exy) * (1 + V1[j][ibxy - 1] * dif + V1[j][ibxy + m1 - 1] * dif1) + (1 - (1 + b) * eb) * (V2[j][ib - 1][0] * dif1 + V2[j][ib + m1 - 1][0] * dif)
                if vn > vmax:
                    vmax = vn
                    bmax = b
            if flag == 1:
                kl = k + l - m2
                b1max = x
                for m in range(k, kl + 1, m1):
                    b = (m - 1) * del_
                    eb = math.exp(-b)
                    exyb = ex * ey / eb
                    ibx = m - k + 1
                    ibxy = k + l - m
                    vn = x * (eb - exy) * (1 + V1[j][ibxy - 1]) + (1 - ex - x * ex) * V2[j][k - 1][ibx - 1]
                    if vn > vmax:
                        b1max = b
                        vmax = vn
                        type_ = 1
                for m in range(1, 20):
                    b = b1max - del1 + m * delf1
                    if b < x:
                        continue
                    eb = math.exp(-b)
                    exyb = ex * ey / eb
                    aux = b / del1
                    ib = int(np.floor(aux))
                    bdif = aux - ib
                    bdif1 = 1 - bdif
                    ib = ib * m1 + 1
                    ibxy = k + l - ib - m1
                    ibx = ib - k + 1
                    vn = x * (eb - exy) * (1 + V1[j][ibxy - 1] * bdif + V1[j][ibxy + m1 - 1] * bdif1) + (1 - ex - x * eb) * (V2[j][k - 1][ibx - 1] * bdif1 + V2[j][k - 1][ibx + m1 - 1] * bdif)
                    if vn > vmax:
                        b1max = b
                        vmax = vn
                        type_ = 1
            if type_ == 0:
                b = bmax
            else:
                b = b1max
            beta[k - 1][l - 1] = b
            V2[i][k - 1][l - 1] = vmax / (1 - ex - x * exy)
            xy = V2[i][k - 1][l - 1] - V2[j][k - 1][l - 1]
            if xy < xmin:
                xmin = xy
            if xy > xmax:
                xmax = xy
            if xmin < 0:
                print("curr V2 = %f, prev V2 = %f, k = %d, l = %d" % (V2[i][k - 1][l - 1], V2[j][k - 1][l - 1], k, l))
                print("Neg xmin: %f" % xmin)
    b = delf1
    for l in range(386, 590):
        x = del_ * (l - 1)
        ex = math.exp(-x)
        vn = 0.0
        vmax = vn
        for b in np.arange(b, x, delf):
            eb = math.exp(-b)
            exb = ex / eb
            aux = b / del1
            ib = int(np.floor(aux))
            bdif = aux - ib
            bdif1 = 1 - bdif
            ib = ib * m1 + 1
            aux = (x - b) / del1
            ibx = int(np.floor(aux))
            bxdif = aux - ibx
            bxdif1 = 1 - bxdif
            ibx = ibx * m1 + 1
            if (ibx < 385):
                vn = (eb - (1 + x - b) * ex) * (V2[j][ibx - 1][0] * bxdif1 + V2[j][ibx + m1 - 1][0] * bxdif) + b * (eb - ex) * (1 + V1[j][ibx - 1] * bxdif1 + V1[j][ibx + m1 - 1] * bxdif) + (1 - (1 + b) * eb) * (V2[j][ib - 1][0] * bdif1 + V2[j][ib + m1 - 1][0] * bdif)
            else:
                vn = (eb - (1 + x - b) * ex) * (V2[j][ibx - 1][0] * bxdif1 + V2[j][ibx + m1 - 1][0] * bxdif) + b * (eb - ex) * (1 + V1[j][ibx - 1]) + (1 - (1 + b) * eb) * (V2[j][ib - 1][0] * bdif1 + V2[j][ib + m1 - 1][0] * bdif)
            if vn > vmax:
                vmax = vn
                bmax = b
        b = bmax
        beta[l - 1][0] = b
        V2[i][l - 1][0] = vmax / (1 - (1 + x) * ex)
        x = V2[i][l - 1][0] - V2[j][l - 1][0]
        if x < xmin:
            xmin = x
        if x > xmax:
            xmax = x
        if xmin < 0:
            print("curr V2 = %f, prev V2 = %f, l = %d" % (V2[i][l - 1][0], V2[j][l - 1][0], l))
            print("ibx=%d, ib =%d" % (ibx, ib))
            print("Neg xmin: %f" % xmin)
out_file2 = open("out_russian_modified.txt", "w")
for i in range(385):
    out_file2.write("%lf, %lf, %lf\n" % (i * del_, V1[1][i], alph[i]))
out_file2.close()
print("%.10lf-%.10lf=%.10lf" % (xmax, xmin, xmax - xmin))

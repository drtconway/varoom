import math
import os
import random
import subprocess
import sys

facCache = {}

def sqr(x):
    return x*x

def factorial(n):
    global facCache
    if n not in facCache:
        f = 1
        j = n
        while j > 1:
            f *= j
            j -= 1
        facCache[n] = f
    return facCache[n]

def pois(lam):
    c = 0.767 - 3.36/lam
    beta = math.pi/math.sqrt(3.0*lam)
    alpha = beta*lam
    k = math.log(c) - lam - math.log(beta)

    while True:
        u = random.random()
        x = (alpha - math.log((1.0 - u)/u))/beta
        n = int(math.floor(x + 0.5))
        if n < 0:
            continue
        v = random.random()
        y = alpha - beta*x
        lhs = y + math.log(v/(1.0 + sqr(math.exp(y))))
        rhs = k + n*math.log(lam) - math.log(factorial(n))
        if lhs <= rhs:
            return n

def readFasta(f):
    nm = None
    seq = []
    for l in f:
        l = l.strip()
        if len(l) and l[0] == '>':
            if nm is not None:
                yield (nm, ''.join(seq))
            nm = l[1:].strip()
            seq = []
        else:
            seq.append(l)
    if nm is not None:
        yield (nm, ''.join(seq))

def openFile(fn):
    if fn == "-":
        return sys.stdin
    if fn.endswith(".gz") and os.path.exists(fn):
        p = subprocess.Popen(['gunzip', '-c', fn],
                             bufsize=1024*1024,
                             stdout=subprocess.PIPE)
        return p.stdout
    return open(fn)

alts = {'A':['C','G','T'],'C':['A','G','T'],'G':['C','A','T'],'T':['C','G','A']}

def mkProbs(E, N, r):
    global alts
    u = random.random()
    s = r
    if u < N:
        s = random.choice(alts[r])
    pr = {}
    if random.random() < 0.5:
        v = 1.0 - E
        pr[r] = v/2.0
        pr[s] = v/2.0
        ss = set(alts[r]) - set([s])
        for b in ss:
            pr[b] = E/2.0
        return pr.items()
    else:
        v = 1.0 - E
        pr[s] = v
        for b in alts[s]:
            pr[b] = E/3.0
        return pr.items()

def mkCov(C, px):
    c = pois(C)
    counts = {'A':0,'C':0,'G':0,'T':0}
    for i in range(c):
        u = random.random()
        for j in range(4):
            if u < px[j][1]:
                counts[px[j][0]] += 1
                break
            u -= px[j][1]
    return counts

bed = []
with open(sys.argv[1]) as f:
    for l in f:
        t = l.split()
        ch = t[0]
        st = int(t[1])
        en = int(t[2])
        bed.append((ch, st, en))

S = int(sys.argv[2])
C = int(sys.argv[3])
E = float(sys.argv[4])
N = float(sys.argv[5])
G = sys.argv[6]

lastCh = None
lastSeq = None

random.seed(S)
for itm in bed:
    (ch, st, en) = itm
    if ch != lastCh:
        lastCh = ch
        with openFile(G + "/" + ch + ".fa.gz") as f:
            for (nm,seq) in readFasta(f):
                lastSeq = seq
    for pos in range(st+1, en+1):
        r = lastSeq[pos].upper()
        pr = mkProbs(E, N, r)
        for (b,c) in sorted(mkCov(C, pr).items()):
            if c > 0:
                print '%s\t%d\t%s\t%d' % (ch, pos, b, c)



def zero_array(n):
    a = []
    for i in range(n):
        a.append(0.0)
    return a


# see: http://blog.plover.com/math/choose.html
def binom(n, k):
    if k == 0:
        return 1
    if n == 0:
        return 0
    m = n
    b = 1
    for j in range(k):
        b = b*m
        b = b/(j + 1)
        m -= 1
    return b


# float to prevent overflow
def fac(n):
    r = 1.0
    for i in range(n):
        r *= float(i + 1)
    return r


# float to prevent overflow
def fac2(n):
    if n < 0:
        r = float(n + 2)
        i = n + 4
        while i < 2:
            r *= float(i)
            i += 2
        return 1.0/r
    elif n == 0:
        return 1.0
    else:
        r = float(n)
        i = n - 2
        while i > 0:
            r *= float(i)
            i -= 2
        return r


def get_cs_trans(max_l_value):
    import math

    n = 0
    for l in range(2, max_l_value + 1):
        t = (l+1)*(l+2)/2
        n += t*t

    cs_trans = []
    cs_trans.append([[1.0]])
    cs_trans.append([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

    for l in range(2, max_l_value + 1):

        nc = (l+1)*(l+2)/2
        ns = 2*l + 1

        legendre_coef = zero_array(l + 1)
        cossin_coef = zero_array((l + 1)*(l + 1))
        tmat = zero_array(int(nc*ns))

        k = 0
        while k <= l/2:
            legendre_coef[l - 2*k] = ((-1.0)**k/2.0**l)*binom(l, k)*binom(2*(l - k), l)
            k += 1

        for m in range(l + 1):
            cossin_coef[m] = 1.0
            for k in range(1, m + 1):
                cossin_coef[k*(l+1) + m] += cossin_coef[(k-1)*(l+1) + m - 1]*((-1.0)**(k-1))
                if m > k:
                    cossin_coef[k*(l+1) + m] += cossin_coef[k*(l+1) + m - 1]

        for m in range(l + 1):
            if m == 0:
                cm = 1.0
            else:
                cm = math.sqrt(2.0*fac(l - m)/fac(l + m))
            cm = cm/math.sqrt(fac2(2*l - 1))

            k = (l - m) % 2
            while k <= (l - m):
                if m > 0:
                    legendre_coef[k] = float(k + 1)*legendre_coef[k+1]
                cmk = cm*legendre_coef[k]
                i = 0
                while i <= (l - k - m)/2:
                    cmki = cmk*binom((l - k - m)/2, i)
                    for j in range(i + 1):
                        cmkij = cmki*binom(i, j)
                        for n in range(m + 1):
                            ix = l - 2*j - m + n
                            ix = ix*(ix + 1)/2 + l + 1 - m - 2*i
                            if n % 2 == 1:
                                ilm = 1 + l - m
                            else:
                                ilm = 1 + l + m
                            tmat[int((ilm - 1)*nc + ix - 1)] += cmkij*cossin_coef[n*(l+1) + m]
                    i += 1
                k += 2

        tc = []
        for i in range(int(nc)):
            ts = []
            for j in range(ns):
                ts.append(tmat[int(j*nc + i)])
            tc.append(ts)
        cs_trans.append(tc)

    return cs_trans

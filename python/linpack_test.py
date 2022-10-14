import numpy as np
import time

def linpack_r_sp(n, niter):
    a = np.random.rand(n, n)
    a = np.array(a, dtype=np.float32)

    b = np.random.rand(n, n)
    b = np.array(b, dtype=np.float32)


    tstart = time.time()
    for i in range(niter):
        c = a @ b
    tend = time.time()

    tot_time = tend - tstart 
    flops = 2.0 * n**3 * niter / (tend-tstart) / 1.0e09

    print("**************************************************")
    print("linpack for real(sp) n = {0:6d}".format(n))
    print("Number of iterations   = {0:6d}".format(niter))
    print("Elasped time in s      = {0:.2f}".format(tot_time))
    print("Perfromance in GFLOPS  = {0:.2f}".format(flops))

    return None


def linpack_r_dp(n, niter):
    a = np.random.rand(n, n)
    a = np.array(a, dtype=np.float64)

    b = np.random.rand(n, n)
    b = np.array(b, dtype=np.float64)


    tstart = time.time()
    for i in range(niter):
        c = a @ b
    tend = time.time()

    tot_time = tend - tstart 
    flops = 2.0 * n**3 * niter / (tend-tstart) / 1.0e09

    print("**************************************************")
    print("linpack for real(dp) n = {0:6d}".format(n))
    print("Number of iterations   = {0:6d}".format(niter))
    print("Elasped time in s      = {0:.2f}".format(tot_time))
    print("Perfromance in GFLOPS  = {0:.2f}".format(flops))

    return None


def linpack_c_sp(n, niter):
    a = np.random.rand(n, n) + 1j*np.random.rand(n, n)
    a = np.array(a, dtype=np.complex64)

    b = np.random.rand(n, n) + 1j*np.random.rand(n, n)
    b = np.array(b, dtype=np.complex64)


    tstart = time.time()
    for i in range(niter):
        c = a @ b
    tend = time.time()

    tot_time = tend - tstart 
    flops = 6.0 * n**3 * niter / (tend-tstart) / 1.0e09

    print("**************************************************")
    print("linpack for complex(sp) n = {0:6d}".format(n))
    print("Number of iterations      = {0:6d}".format(niter))
    print("Elasped time in s         = {0:.2f}".format(tot_time))
    print("Perfromance in GFLOPS     = {0:.2f}".format(flops))

    return None


def linpack_c_dp(n, niter):
    a = np.random.rand(n, n) + 1j*np.random.rand(n, n)
    a = np.array(a, dtype=np.complex128)

    b = np.random.rand(n, n) + 1j*np.random.rand(n, n)
    b = np.array(b, dtype=np.complex128)


    tstart = time.time()
    for i in range(niter):
        c = a @ b
    tend = time.time()

    tot_time = tend - tstart 
    flops = 6.0 * n**3 * niter / (tend-tstart) / 1.0e09

    print("**************************************************")
    print("linpack for complex(dp) n = {0:6d}".format(n))
    print("Number of iterations      = {0:6d}".format(niter))
    print("Elasped time in s         = {0:.2f}".format(tot_time))
    print("Perfromance in GFLOPS     = {0:.2f}".format(flops))

    return None


n = 500
niter = 200
linpack_r_sp(n, niter)
linpack_r_dp(n, niter)
linpack_c_sp(n, niter)
linpack_c_dp(n, niter)
print("**************************************************")
print()
print()
print()

n = 1000
niter = 100
linpack_r_sp(n, niter)
linpack_r_dp(n, niter)
linpack_c_sp(n, niter)
linpack_c_dp(n, niter)
print("**************************************************")
print()
print()
print()

n = 2000
niter = 10
linpack_r_sp(n, niter)
linpack_r_dp(n, niter)
linpack_c_sp(n, niter)
linpack_c_dp(n, niter)
print("**************************************************")
print()
print()
print()

n = 3500
niter = 5
linpack_r_sp(n, niter)
linpack_r_dp(n, niter)
linpack_c_sp(n, niter)
linpack_c_dp(n, niter)
print("**************************************************")
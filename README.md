# liblsoda

This is a C library that implements the LSODA algorithm from Linda Petzold and Alan Hindmarsh, which solves the initial value problem for stiff or nonstiff systems of first order ordinary differential equations

```
dy/dt = f(t,y)
```

or, in component form,

```
dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
```

## Building

This code has only been tested on Linux to date.

```
git clone http://github.com/sdwfrost/liblsoda
cd liblsoda
make
```

If you want to compile with OpenMP, using `CFLAGS=-fopenmp make` instead.

## Testing

```
cd test
make
export LD_LIBRARY_PATH=../src:$LD_LIBRARY_PATH
./test
```

This should give the following.

```
at t=   4.0000e-01 y=   9.851712e-01   3.386380e-05   1.479493e-02
at t=   4.0000e+00 y=   9.055333e-01   2.240655e-05   9.444430e-02
at t=   4.0000e+01 y=   7.158403e-01   9.186334e-06   2.841505e-01
at t=   4.0000e+02 y=   4.505250e-01   3.222964e-06   5.494717e-01
at t=   4.0000e+03 y=   1.831976e-01   8.941773e-07   8.168015e-01
at t=   4.0000e+04 y=   3.898729e-02   1.621940e-07   9.610125e-01
at t=   4.0000e+05 y=   4.936362e-03   1.984221e-08   9.950636e-01
at t=   4.0000e+06 y=   5.161833e-04   2.065787e-09   9.994838e-01
at t=   4.0000e+07 y=   5.179804e-05   2.072027e-10   9.999482e-01
at t=   4.0000e+08 y=   5.283676e-06   2.113481e-11   9.999947e-01
at t=   4.0000e+09 y=   4.658665e-07   1.863467e-12   9.999995e-01
at t=   4.0000e+10 y=   1.431150e-08   5.724606e-14   1.000000e+00
```


## Contributors

This code has a compllicated history; all I have done is to refactor the code so that it compiles cleanly, and generates a shared library.

- The original code by Petzold and Marsh was written in Fortran, then converted to C. It is available from [here](http://www.ccl.net/cca/software/SOURCES/C/kinetics2/index.shtml), and the code is mirrored in the `orig/lsoda-orig` subdirectory.
- Heng Li [@lh3](http://github.com/lh3) refactored the code, which is available [here](http://github.com/lh3/misc). The code is mirrored in the `orig/lsoda-lh3` subdirectory.
- Yu Feng [@rainwoodman](http://github.com/rainwoodman) refactored Heng Li's code to make it reentrant. The code is available [here](https://github.com/rainwoodman/psphray), and is mirrored in the `orig/lsoda-rainwoodman` subdirectory.

## References

1.  Alan C. Hindmarsh.  ODEPACK, a systematized collection of ODE solvers. In: Scientific Computing, R. S. Stepleman et al. (eds.), North-Holland, Amsterdam, 1983, pp. 55-64. A preprint is available [here](https://computation.llnl.gov/casc/nsde/pubs/u88007.pdf).
2.  Linda R. Petzold. Automatic selection of methods for solving stiff and nonstiff systems of ordinary differential equations. SIAM J. Sci. Stat. Comput. 4 (1983), pp. 136-148.

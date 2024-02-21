# Submodule Integrate

There are implemented several routines for integration of functions, and also sampled values via Trapezoid and Simpson methods. 

Each algorithm computes an approximation to a definite integral of the
form,

\f[J \approx I := \int_a^b f(x) w(x) dx\f]

where \f$w(x)\f$ is a weight function, which may be \f$w(x) = 1\f$.


## @anchor integwhichfunction Integrable functions

All routines described here may integrate functions that return either **real** or **complex** values. They may have one of the following signatures:

* real function `f(x)` of a real argument `x`
* complex function `f(x)` of a real argument `x`
* real function `f(x,args)` of real arguments `x`, and array `args`
* complex function `f(x,args)` of real arguments `x`, and array `args`

For functions that need extra arguments besides the integration variable `x`, the argument array must be passed to the integration routines.

There is also the possibility of using
- [type(nf_rfunction)](@ref func_integ::nf_rfunction) for real functions or
- [type(nf_cfunction)](@ref func_integ::nf_cfunction) for complex functions



## @anchor desctrapezoid Trapezoid and Simpson rules


There are implemented routines for integration of functions by simple [Trapezoid and Simpson](@ref qsimpson) rules via two routines [trapz](@ref qsimpson::trapz) and [simps](@ref qsimpson::simps). The use of both routines is very similar. For instance:

If we have data sampled to a grid of **ordered** values, we can integrate them by using trapezoid or simpson quadrature similarly as with functions. In order to use them we have only to call [trapz](@ref qsimpson::trapz) or [simps](@ref qsimpson::simps)

@include ex_integr_trapz.f90

The use of Simpson routine is similar, and need only to replace `trapz()` by  `simps()`

## @anchor desciadsimps Adaptive Simpson method ##

In addition to The simple, fixed points, Simpson quadrature routines there are implemented Globally Adaptive Simpson routines for [finite](@ref qadaptive::iads) and [infinite](@ref qadaptive::iads) integration domains.

Routine [iads](@ref qadaptive::iads) follows a scheme in which the integration interval is bisected, and the Simpson integration result for the sum of the subintervals is compared with the result for the parent interval. For each subdivision, only three additional function evaluations are needed. This routine is regarded as a good choice when moderate precision (\f$10^{-3}, 10^{-4}\f$) is required.

@include ex_integr_iads.f90


Routine [iadsi](@ref qadaptive::iadsi) is available for semi-infinite integration domains.

For instance, for the function

@image html ex_integr_iadsi.png

we see that the main contribution comes from the region below \f$x \sim 1\f$. Then, we can use that information as argument `brkpts`, as shown in the following example

@include ex_integr_iadsi.f90

## @anchor descthsh Integration using the Tanh-sinh scheme ##

Routine [qnthsh](@ref qtanhsinh::qnthsh) is an implementation of the Tanh-sinh quadrature method.
As explained in [Wikipedia](https://en.wikipedia.org/wiki/Tanh-sinh_quadrature), the *tanh-sinh* or *Double exponential* quadrature is a method for numerical integration introduced by [Hidetosi Takahasi and Masatake Mori in 1974](https://doi.org/10.2977%2Fprims%2F1195192451). It is based in the change of variables

\f[ x := g(t) = \tanh\left(\frac{\pi}{2}\sinh t\right)\,\f]

to transform an integral on the interval x ∈ (−1, +1) to an integral on the entire real line t ∈ (−∞, +∞). After this transformation, the integrand decays with a double exponential rate, and thus, this method is also known as the "Double Exponential (DE) formula" ([Mori, Masatake (2005), "Discovery of the Double Exponential Transformation and Its Developments", Publications of the Research Institute for Mathematical Sciences, 41 (4): 897-935](https://doi.org/10.2977/prims/1145474600)).

For a given step size h, the integral is approximated by the sum
\f[ \int_{-1}^{1} f(x) \, dx \approx \sum_{k=-\infty}^\infty w_{k} f(x_{k}),\f]

with the abscissas
\f[ x_k := g(h k) = \tanh\left(\frac{\pi}{2}\sinh kh\right) \f]

and the weights
\f[w_k := g'( h k) = \frac{\frac{\pi}{2}h \cosh(kh)}{\cosh^{2} \big(\frac{\pi}{2} \sinh (kh) \big)}.\f]

The Tanh-Sinh method is quite insensitive to endpoint behavior. Should singularities or infinite derivatives exist at one or both endpoints of the (−1, +1) interval, these are mapped to the (−∞,+∞) endpoints of the transformed interval, forcing the endpoint singularities and infinite derivatives to vanish. This results in a great enhancement of the accuracy of the numerical integration procedure, which is typically performed by the Trapezoidal Rule. In most cases, the transformed integrand displays a rapid roll-off (decay), enabling the numerical integrator to quickly achieve convergence. 

### Examples of use ###

The most simple example of using the routine is given by

@include{lineno} ex_integr_thsh.f90

Here in the call in line 10 is the simplest possible and only uses the required arguments.

## @anchor descquadpack Quadpack routines ##

The code implemented here has been slightly adapted from the original **[Quadpack](http://www.netlib.org/quadpack/)** routines (by Piessens, de Doncker-Kapenga, Ueberhuber and Kahaner). It has been ported to Modern Fortran by following several sources: [Gnu Scientific Library](https://www.gnu.org/software/gsl/), [J. Burkardt's page](https://people.sc.fsu.edu/~jburkardt/f_src/quadpack_double/quadpack_double.html), and [SciFortran](https://github.com/aamaricci/SciFortran). Moreover, the routines have been tweaked to work with both real and complex functions.

There are routines for adaptive and non-adaptive integration of general functions, with specific routines for some particular difficult cases. These include integration over infinite and semi-infinite ranges, singular integrals, including potential and logarithmic singularities, computation of Cauchy principal values and oscillatory integrals. 

The algorithms in *Quadpack* use a naming convention based on the following letters:

| Symbol | Explanation                                  |
|--------|----------------------------------------------|
| Q      | quadrature routine                           |
| N      | non-adaptive integrator                      |
| A      | adaptive integrator                          |
| G      | general integrand (user-defined)             |
| W      | weight function with integrand               |
| S      | singularities can be more readily integrated |
| P      | points of special difficulty can be supplied |
| I      | infinite range of integration                |
| O      | oscillatory weight function, cos or sin      |
| F      | Fourier integral                             |
| C      | Cauchy principal value                       |

The algorithms are built on pairs of quadrature rules, a higher order
rule and a lower order rule.  The higher order rule is used to compute
the best approximation to an integral over a small range.  The
difference between the results of the higher order rule and the lower
order rule gives an estimate of the error in the approximation.

All subroutines present a similar interface. They accept several optional input and output arguments but, in their simpler form, all integrators without weight functions may be called as

```{.f90}
  call integrator(f, a, b, IntVal)
```
where `f` is a real or complex function `f(x,[args])` as described in [Integrable functions](@ref integwhichfunction), `IntVal` is the value returned with an approximation to the integral of `f` over the (possibly infinite) interval `(a,b)`.

There are available the following integration routines:

### Subroutine QNG for non-adaptive integration ###

The [subroutine qng](@ref quadpack::qng) is a simple non-adaptive automatic integrator, based on a sequence of rules with increasing degree of algebraic precision (Patterson, 1968). It applies the Gauss-Kronrod 10-point, 21-point, 43-point and 87-point integration rules in succession until an estimate of the integral is achieved within the desired absolute and relative error limits. Each successive approximation reuses the previous function evaluations.

#### Examples of use ####

@include{lineno} ex_integr_qng.f90

### Subroutine QAG for simple adaptive integration ###

The subroutine [qag](@ref quadpack::qag) uses a global adaptive scheme in which the interval is divided and at each step the integral is calculated for each subinterval, and the one with the larger error is bisected. This scheme is recusively applied and the most difficult points are calculated in finer grids.

#### Examples of use ####

@include{lineno} ex_integr_qag.f90


### Subroutine QAGS for general integration with singularities ###

Subroutine [qags](@ref quadpack::qags) is a global adaptive algorithm with bisection of the integration interval and a Wynn-epsilon algorithm to extrapolate and speed-up the convergence of the integral. This routine is able to integrate functions with singularities of many types. 

The limits of integration may be finite or infinite. By default it uses a Gauss-Kronrod rule of 21 points for finite limits, and a 15-points rule for infinite limits. The rule may be given to the routine as an argument.

The integration over the infinite intervals are is performed by means of mappings to the semi-open interval \f$(0,1]\f$.

#### Examples of use ####

@include{lineno} ex_integr_qags.f90

### Subroutine QAGP integration with known singular points ###

Subroutine [qagp](@ref quadpack::qagp) uses the same strategy than [qags](@ref quadpack::qags) to integrate a function that presents singularities in the interior of the integration domain. The user has to supply an array with the points where difficulties arise. 

#### Examples of use ####

@include{lineno} ex_integr_qagp.f90

### Subroutine QAWO for Oscillatory functions ###

Subroutine [qawo](@ref quadpack::qawo) computes the integral of a function `f` multiplied by a oscillatory function of \f$\cos{(\omega x)}\f$ (for `flgw`=1) or \f$\sin{(\omega x)}\f$ (for `flgw`=1), 

#### Examples of use ####

@include{lineno} ex_integr_qawo.f90

### Subroutine QAWF for Fourier Transforms ###

Subroutine [qawf](@ref quadpack::qawf) computes the (possibly semi-infinite) cosine (for `flgw`=1) or sine (for `flgw`=2) Fourier transform of `f` using an adaptive algorithm

#### Examples of use ####

@include{lineno} ex_integr_qawf.f90

### Subroutine QAWS for integration of singular functions ###

Subroutine [qaws](@ref quadpack::qaws) computes the integral of a function `f` multiplied by a weight function with an algebraic-logarithmic singularity chosen according to the flag `flgw`,

| `flgw` | Weight function                                            |
|--------|------------------------------------------------------------|
| 1      | \f$(x - a)^\alpha (b - x)^\beta\f$                         |
| 2      | \f$(x - a)^\alpha (b - x)^\beta \log{(x-a)}\f$             |
| 3      | \f$(x - a)^\alpha (b - x)^\beta \log{(b-x)}\f$             |
| 4      | \f$(x - a)^\alpha (b - x)^\beta \log{(x-a)} \log{(b-x)}\f$ |

The points \f$(a,b)\f$ are the limits of the integration domain.

#### Examples of use ####

@include{lineno} ex_integr_qaws.f90

### Subroutine QAWC for Cauchy principal values ###

Subroutine [qawc](@ref quadpack::qawc) computes Cauchy Principal Value of the integral of `f` using the bisection strategy of [qag](@ref quadpack::qag) 

#### Examples of use ####

@include{lineno} ex_integr_qawc.f90

### @anchor qp_info Use of the extra parameter `info`

Most adaptive routines allow for an optional input/output argument `info`. It is mainly used to pass information on details of the calculations. This argument is a variable of the type [d_qp_extra](@ref quadpack::d_qp_extra) for integration of real functions, and of the type [c_qp_extra](@ref quadpack::c_qp_extra) for integration of complex functions.

Following Piessend *et al* we illustrate one of its use considering the integral 
\f[
\int_{0}^{1} \left| x^{2} + 2 x -1 \right|^{-1/2} \, dx = \arcsin(1/\sqrt(3)) - \frac{\pi}{2} - \log(3)/2 \approx 1.504622
\f]

In this case the integral converges

@include{lineno} ex_integr_info.f90

giving as result information on the algorithm progress

```{shell-session}
integrate(1/sqrt(|x^2 + 2x -2|), x, 0, 1) = 1.5045599601491 (Error=0.0001175019052) with N = 735
Difference = 6.2802309494403e-05 relative to analytical value
Error code = 0
Meaning:
-------
Routine qags:  Normal and reliable termination of the routine.  It is assumed that the requested accuracy has been achieved.

 j        Limits of subinterval S(j)     Integral over S(j)    Error on S(j)    Error relative 
 1 7.32050657272E-01  7.32050895691E-01  6.61788951586E-04  2.67199722355E-04  4.03753676628E-01
 2 7.50000000000E-01  8.75000000000E-01  2.59695449021E-01  2.34942884257E-05  9.04686181996E-05
 3 7.32055664062E-01  7.32116699219E-01  6.35457321524E-03  2.21241348822E-05  3.48160830520E-03
 4 7.34375000000E-01  7.50000000000E-01  9.20419410896E-02  6.43439223180E-06  6.99071766157E-05
 5 7.32050895691E-01  7.32051849365E-01  7.77806648660E-04  1.18010670163E-06  1.51722372605E-03
 6 7.32421875000E-01  7.34375000000E-01  3.10999052060E-02  3.36607750891E-07  1.08234333405E-05
 7 7.31445312500E-01  7.31933593750E-01  1.48085753181E-02  2.20888795975E-08  1.49162759570E-06
 8 0.00000000000E+00  5.00000000000E-01  4.31717842526E-01  7.65314094754E-10  1.77271824179E-09
 9 7.32051849365E-01  7.32055664062E-01  1.27128139071E-03  5.82935274600E-10  4.58541499043E-07
10 6.87500000000E-01  7.18750000000E-01  1.03290479678E-01  5.00387915491E-10  4.84447276312E-09
11 7.26562500000E-01  7.30468750000E-01  3.68841542687E-02  3.09541720462E-10  8.39226834936E-09
12 7.32040405273E-01  7.32048034668E-01  1.67639077543E-03  4.48780642525E-11  2.67706461466E-08
13 7.32177734375E-01  7.32421875000E-01  8.59296586961E-03  4.15971773602E-12  4.84084051903E-10
14 7.32048034668E-01  7.32049942017E-01  7.89650429625E-04  1.86687676455E-12  2.36418128137E-09
15 6.25000000000E-01  6.87500000000E-01  1.26121798856E-01  1.15514630686E-12  9.15897424031E-12
16 7.30468750000E-01  7.31445312500E-01  1.63018589630E-02  9.07839605128E-13  5.56893301058E-11
17 7.18750000000E-01  7.26562500000E-01  4.43801352397E-02  4.94821770138E-13  1.11496228541E-11
18 5.00000000000E-01  6.25000000000E-01  1.70177840734E-01  1.29207244330E-13  7.59248347332E-13
19 7.32025146484E-01  7.32040405273E-01  1.97766297002E-03  3.28124217391E-14  1.65915134361E-11
20 7.32050418854E-01  7.32050657272E-01  2.53372044356E-04  1.18649893348E-14  4.68283285354E-11
21 7.31994628906E-01  7.32025146484E-01  2.61075315192E-03  2.57842783013E-15  9.87618392123E-13
22 8.75000000000E-01  1.00000000000E+00  1.45769659141E-01  1.76224733434E-15  1.20892601707E-14
23 7.31933593750E-01  7.31994628906E-01  3.57974836972E-03  1.00512858960E-15  2.80781911405E-13
24 7.32049942017E-01  7.32050418854E-01  3.29764477394E-04  3.75908490859E-16  1.13993021271E-12
25 7.32116699219E-01  7.32177734375E-01  3.38357303442E-03  9.07748753959E-17  2.68281117247E-14
```

As we can see, the general integrator `qags` is able to integrate this singular function but needs several function evaluations, with an important error in the vicinity of \f$\sqrt{3}-1\f$. We can use this additional information and replace the routine `qags` by `qagp`. It only need to change line number XX to:
```{.f90}
call qagp(fquad4510, Zero, 1._dp, [sqrt(3._dp) - 1], Integ1, epsabs=Zero, epsrel=1.e-4_dp, abserr=err, Neval=N, ier=ier, info=info)
```

obtaining a similar (and somewhat better) result with less subdivisions of the interval and function evaluations

```{shell-session}
integrate(1/sqrt(|x^2 + 2x -2|), x, 0, 1) = 1.5046227555966 (Error=3.5719196844752e-07) with N = 462
Difference = 6.8619860904562e-09 relative to analytical value
Error code = 0
Meaning:
-------
Routine qagp: Normal and reliable termination of the routine.  It is assumed that the requested accuracy has been achieved.

 j        Limits of subinterval S(j)     Integral over S(j)    Error on S(j)    Error relative 
 1 7.09174219832E-01  7.32050807569E-01  1.60069187812E-01  7.73744297312E-02  4.83381160290E-01
 2 7.32050807569E-01  7.40424219832E-01  9.66937903941E-02  4.68487167016E-02  4.84505949251E-01
 3 0.00000000000E+00  3.66025403784E-01  2.93171381140E-01  8.17539976481E-09  2.78860771915E-08
 4 3.66025403784E-01  5.49038105677E-01  1.98297261193E-01  5.21098430467E-09  2.62786499083E-08
 5 8.66025403784E-01  1.00000000000E+00  1.58478399902E-01  3.86710284589E-09  2.44014505969E-08
 6 5.49038105677E-01  6.40544456623E-01  1.37342006395E-01  3.52212881075E-09  2.56449494455E-08
 7 7.99038105677E-01  8.66025403784E-01  1.13597226144E-01  2.80948755739E-09  2.47320084544E-08
 8 6.40544456623E-01  6.86297632096E-01  9.61465434121E-02  2.43826866463E-09  2.53599201604E-08
 9 7.65544456623E-01  7.99038105677E-01  8.08861095193E-02  2.01475142213E-09  2.49084970721E-08
10 6.86297632096E-01  7.09174219832E-01  6.76512094521E-02  1.70646046002E-09  2.52243895393E-08
11 7.48797632096E-01  7.65544456623E-01  5.73966267734E-02  1.43490327440E-09  2.49997840476E-08
12 7.40424219832E-01  7.48797632096E-01  4.06573774170E-02  1.01831363346E-09  2.50462203455E-08
```





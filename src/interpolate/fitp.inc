subroutine percur(iopt, m, x, y, w, k, s, nest, n, t, c, fp,&
  &wrk, lwrk, iwrk, ier)
  !  given the set of data points (x(i),y(i)) and the set of positive
  !  numbers w(i),i=1,2,...,m-1, subroutine percur determines a smooth
  !  periodic spline approximation of degree k with period per=x(m)-x(1).
  !  if iopt=-1 percur calculates the weighted least-squares periodic
  !  spline according to a given set of knots.
  !  if iopt>=0 the number of knots of the spline s(x) and the position
  !  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
  !  ness of s(x) is then achieved by minimalizing the discontinuity
  !  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
  !  n-k-1. the amount of smoothness is determined by the condition that
  !  f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
  !  negative constant, called the smoothing factor.
  !  the fit s(x) is given in the b-spline representation (b-spline coef-
  !  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
  !  subroutine splev.
  !
  !  calling sequence:
  !     call percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,wrk,
  !    * lwrk,iwrk,ier)
  !
  !  parameters:
  !   iopt  : integer flag. on entry iopt must specify whether a weighted
  !           least-squares spline (iopt=-1) or a smoothing spline (iopt=
  !           0 or 1) must be determined. if iopt=0 the routine will start
  !           with an initial set of knots t(i)=x(1)+(x(m)-x(1))*(i-k-1),
  !           i=1,2,...,2*k+2. if iopt=1 the routine will continue with
  !           the knots found at the last call of the routine.
  !           attention: a call with iopt=1 must always be immediately
  !           preceded by another call with iopt=1 or iopt=0.
  !           unchanged on exit.
  !   m     : integer. on entry m must specify the number of data points.
  !           m > 1. unchanged on exit.
  !   x     : real array of dimension at least (m). before entry, x(i)
  !           must be set to the i-th value of the independent variable x,
  !           for i=1,2,...,m. these values must be supplied in strictly
  !           ascending order. x(m) only indicates the length of the
  !           period of the spline, i.e per=x(m)-x(1).
  !           unchanged on exit.
  !   y     : real array of dimension at least (m). before entry, y(i)
  !           must be set to the i-th value of the dependent variable y,
  !           for i=1,2,...,m-1. the element y(m) is not used.
  !           unchanged on exit.
  !   w     : real array of dimension at least (m). before entry, w(i)
  !           must be set to the i-th value in the set of weights. the
  !           w(i) must be strictly positive. w(m) is not used.
  !           see also further comments. unchanged on exit.
  !   k     : integer. on entry k must specify the degree of the spline.
  !           1<=k<=5. it is recommended to use cubic splines (k=3).
  !           the user is strongly dissuaded from choosing k even,together
  !           with a small s-value. unchanged on exit.
  !   s     : real.on entry (in case iopt>=0) s must specify the smoothing
  !           factor. s >=0. unchanged on exit.
  !           for advice on the choice of s see further comments.
  !   nest  : integer. on entry nest must contain an over-estimate of the
  !           total number of knots of the spline returned, to indicate
  !           the storage space available to the routine. nest >=2*k+2.
  !           in most practical situation nest=m/2 will be sufficient.
  !           always large enough is nest=m+2*k,the number of knots needed
  !           for interpolation (s=0). unchanged on exit.
  !   n     : integer.
  !           unless ier = 10 (in case iopt >=0), n will contain the
  !           total number of knots of the spline approximation returned.
  !           if the computation mode iopt=1 is used this value of n
  !           should be left unchanged between subsequent calls.
  !           in case iopt=-1, the value of n must be specified on entry.
  !   t     : real array of dimension at least (nest).
  !           on successful exit, this array will contain the knots of the
  !           spline,i.e. the position of the interior knots t(k+2),t(k+3)
  !           ...,t(n-k-1) as well as the position of the additional knots
  !           t(1),t(2),...,t(k+1)=x(1) and t(n-k)=x(m),..,t(n) needed for
  !           the b-spline representation.
  !           if the computation mode iopt=1 is used, the values of t(1),
  !           t(2),...,t(n) should be left unchanged between subsequent
  !           calls. if the computation mode iopt=-1 is used, the values
  !           t(k+2),...,t(n-k-1) must be supplied by the user, before
  !           entry. see also the restrictions (ier=10).
  !   c     : real array of dimension at least (nest).
  !           on successful exit, this array will contain the coefficients
  !           c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x)
  !   fp    : real. unless ier = 10, fp contains the weighted sum of
  !           squared residuals of the spline approximation returned.
  !   wrk   : real array of dimension at least (m*(k+1)+nest*(8+5*k)).
  !           used as working space. if the computation mode iopt=1 is
  !           used, the values wrk(1),...,wrk(n) should be left unchanged
  !           between subsequent calls.
  !   lwrk  : integer. on entry,lwrk must specify the actual dimension of
  !           the array wrk as declared in the calling (sub)program. lwrk
  !           must not be too small (see wrk). unchanged on exit.
  !   iwrk  : integer array of dimension at least (nest).
  !           used as working space. if the computation mode iopt=1 is
  !           used,the values iwrk(1),...,iwrk(n) should be left unchanged
  !           between subsequent calls.
  !   ier   : integer. unless the routine detects an error, ier contains a
  !           non-positive value on exit, i.e.
  !    ier=0  : normal return. the spline returned has a residual sum of
  !             squares fp such that abs(fp-s)/s <= tol with tol a relat-
  !             ive tolerance set to 0.001 by the program.
  !    ier=-1 : normal return. the spline returned is an interpolating
  !             periodic spline (fp=0).
  !    ier=-2 : normal return. the spline returned is the weighted least-
  !             squares constant. in this extreme case fp gives the upper
  !             bound fp0 for the smoothing factor s.
  !    ier=1  : error. the required storage space exceeds the available
  !             storage space, as specified by the parameter nest.
  !             probably causes : nest too small. if nest is already
  !             large (say nest > m/2), it may also indicate that s is
  !             too small
  !             the approximation returned is the least-squares periodic
  !             spline according to the knots t(1),t(2),...,t(n). (n=nest)
  !             the parameter fp gives the corresponding weighted sum of
  !             squared residuals (fp>s).
  !    ier=2  : error. a theoretically impossible result was found during
  !             the iteration process for finding a smoothing spline with
  !             fp = s. probably causes : s too small.
  !             there is an approximation returned but the corresponding
  !             weighted sum of squared residuals does not satisfy the
  !             condition abs(fp-s)/s < tol.
  !    ier=3  : error. the maximal number of iterations maxit (set to 20
  !             by the program) allowed for finding a smoothing spline
  !             with fp=s has been reached. probably causes : s too small
  !             there is an approximation returned but the corresponding
  !             weighted sum of squared residuals does not satisfy the
  !             condition abs(fp-s)/s < tol.
  !    ier=10 : error. on entry, the input data are controlled on validity
  !             the following restrictions must be satisfied.
  !             -1<=iopt<=1, 1<=k<=5, m>1, nest>2*k+2, w(i)>0,i=1,...,m-1
  !             x(1)<x(2)<...<x(m), lwrk>=(k+1)*m+nest*(8+5*k)
  !             if iopt=-1: 2*k+2<=n<=min(nest,m+2*k)
  !                         x(1)<t(k+2)<t(k+3)<...<t(n-k-1)<x(m)
  !                       the schoenberg-whitney conditions, i.e. there
  !                       must be a subset of data points xx(j) with
  !                       xx(j) = x(i) or x(i)+(x(m)-x(1)) such that
  !                         t(j) < xx(j) < t(j+k+1), j=k+1,...,n-k-1
  !             if iopt>=0: s>=0
  !                         if s=0 : nest >= m+2*k
  !             if one of these conditions is found to be violated,control
  !             is immediately repassed to the calling program. in that
  !             case there is no approximation returned.
  !
  !  further comments:
  !   by means of the parameter s, the user can control the tradeoff
  !   between closeness of fit and smoothness of fit of the approximation.
  !   if s is too large, the spline will be too smooth and signal will be
  !   lost ; if s is too small the spline will pick up too much noise. in
  !   the extreme cases the program will return an interpolating periodic
  !   spline if s=0 and the weighted least-squares constant if s is very
  !   large. between these extremes, a properly chosen s will result in
  !   a good compromise between closeness of fit and smoothness of fit.
  !   to decide whether an approximation, corresponding to a certain s is
  !   satisfactory the user is highly recommended to inspect the fits
  !   graphically.
  !   recommended values for s depend on the weights w(i). if these are
  !   taken as 1/d(i) with d(i) an estimate of the standard deviation of
  !   y(i), a good s-value should be found in the range (m-sqrt(2*m),m+
  !   sqrt(2*m)). if nothing is known about the statistical error in y(i)
  !   each w(i) can be set equal to one and s determined by trial and
  !   error, taking account of the comments above. the best is then to
  !   start with a very large value of s ( to determine the least-squares
  !   constant and the corresponding upper bound fp0 for s) and then to
  !   progressively decrease the value of s ( say by a factor 10 in the
  !   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
  !   approximation shows more detail) to obtain closer fits.
  !   to economize the search for a good s-value the program provides with
  !   different modes of computation. at the first call of the routine, or
  !   whenever he wants to restart with the initial set of knots the user
  !   must set iopt=0.
  !   if iopt=1 the program will continue with the set of knots found at
  !   the last call of the routine. this will save a lot of computation
  !   time if percur is called repeatedly for different values of s.
  !   the number of knots of the spline returned and their location will
  !   depend on the value of s and on the complexity of the shape of the
  !   function underlying the data. but, if the computation mode iopt=1
  !   is used, the knots returned may also depend on the s-values at
  !   previous calls (if these were smaller). therefore, if after a number
  !   of trials with different s-values and iopt=1, the user can finally
  !   accept a fit as satisfactory, it may be worthwhile for him to call
  !   percur once more with the selected value for s but now with iopt=0.
  !   indeed, percur may then return an approximation of the same quality
  !   of fit but with fewer knots and therefore better if data reduction
  !   is also an important objective for the user.
  !
  !  other subroutines required:
  !    fpbacp,fpbspl,fpchep,fpperi,fpdisc,fpgivs,fpknot,fprati,fprota
  !
  !  references:
  !   dierckx p. : algorithms for smoothing data with periodic and
  !                parametric splines, computer graphics and image
  !                processing 20 (1982) 171-184.
  !   dierckx p. : algorithms for smoothing data with periodic and param-
  !                etric splines, report tw55, dept. computer science,
  !                k.u.leuven, 1981.
  !   dierckx p. : curve and surface fitting with splines, monographs on
  !                numerical analysis, oxford university press, 1993.
  !
  !  author:
  !    p.dierckx
  !    dept. computer science, k.u. leuven
  !    celestijnenlaan 200a, b-3001 heverlee, belgium.
  !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
  !
  !  creation date : may 1979
  !  latest update : march 1987
  !
  !  ..
  !  ..scalar arguments..
  real(8) :: s, fp
  integer :: iopt, m, k, nest, n, lwrk, ier
  !  ..array arguments..
  real(8) :: x(m), y(m), w(m), t(nest), c(nest), wrk(lwrk)
  integer :: iwrk(nest)
  !  ..local scalars..
  real(8) :: per, tol
  integer :: i, ia1, ia2, ib, ifp, ig1, ig2, iq, iz, i1, i2, j1, j2, k1, k2, lwest,&
    &maxit, m1, nmin
  !  ..subroutine references..
  !    perper,pcheck
  !  ..
  !  we set up the parameters tol and maxit
  maxit = 20
  tol = 0.1e-02
  !  before starting computations a data check is made. if the input data
  !  are invalid, control is immediately repassed to the calling program.
  ier = 10
  if (k <= 0 .or. k > 5) return
  k1 = k + 1
  k2 = k1 + 1
  if (iopt < (-1) .or. iopt > 1) return
  nmin = 2 * k1
  if (m < 2 .or. nest < nmin) return
  lwest = m * k1 + nest * (8 + 5 * k)
  if (lwrk < lwest) return
  m1 = m - 1
  do i = 1, m1
    if (x(i) >= x(i + 1) .or. w(i) <= 0.) return
  end do

  if (iopt >= 0) go to 30
  if (n <= nmin .or. n > nest) return
  per = x(m) - x(1)
  j1 = k1
  t(j1) = x(1)
  i1 = n - k
  t(i1) = x(m)
  j2 = j1
  i2 = i1
  do i = 1, k
    i1 = i1 + 1
    i2 = i2 - 1
    j1 = j1 + 1
    j2 = j2 - 1
    t(j2) = t(i2) - per
    t(i1) = t(j1) + per
  end do

  call fpchep(x, m, t, n, k, ier)
  if (ier == 0) go to 40
  return
30 if (s < 0.) return
  if (s == 0. .and. nest < (m + 2 * k)) return
  ier = 0
  ! we partition the working space and determine the spline approximation.
40 ifp = 1
  iz = ifp + nest
  ia1 = iz + nest
  ia2 = ia1 + nest * k1
  ib = ia2 + nest * k
  ig1 = ib + nest * k2
  ig2 = ig1 + nest * k2
  iq = ig2 + nest * k1
  call fpperi(iopt, x, y, w, m, k, s, nest, tol, maxit, k1, k2, n, t, c, fp,&
    &wrk(ifp), wrk(iz), wrk(ia1), wrk(ia2), wrk(ib), wrk(ig1), wrk(ig2),&
    &wrk(iq), iwrk, ier)
end subroutine percur

subroutine curfit(iopt, m, x, y, w, xb, xe, k, s, nest, n, t, c, fp, wrk, lwrk, iwrk, ier)
  !  given the set of data points (x(i),y(i)) and the set of positive
  !  numbers w(i),i=1,2,...,m,subroutine curfit determines a smooth spline
  !  approximation of degree k on the interval xb <= x <= xe.
  !  if iopt=-1 curfit calculates the weighted least-squares spline
  !  according to a given set of knots.
  !  if iopt>=0 the number of knots of the spline s(x) and the position
  !  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
  !  ness of s(x) is then achieved by minimalizing the discontinuity
  !  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
  !  n-k-1. the amount of smoothness is determined by the condition that
  !  f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
  !  negative constant, called the smoothing factor.
  !  the fit s(x) is given in the b-spline representation (b-spline coef-
  !  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
  !  subroutine splev.
  !
  !  calling sequence:
  !     call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
  !
  !  parameters:
  !   iopt  : integer flag. on entry iopt must specify whether a weighted
  !           least-squares spline (iopt=-1) or a smoothing spline (iopt=
  !           0 or 1) must be determined. if iopt=0 the routine will start
  !           with an initial set of knots t(i)=xb, t(i+k+1)=xe, i=1,2,...
  !           k+1. if iopt=1 the routine will continue with the knots
  !           found at the last call of the routine.
  !           attention: a call with iopt=1 must always be immediately
  !           preceded by another call with iopt=1 or iopt=0.
  !           unchanged on exit.
  !   m     : integer. on entry m must specify the number of data points.
  !           m > k. unchanged on exit.
  !   x     : real array of dimension at least (m). before entry, x(i)
  !           must be set to the i-th value of the independent variable x,
  !           for i=1,2,...,m. these values must be supplied in strictly
  !           ascending order. unchanged on exit.
  !   y     : real array of dimension at least (m). before entry, y(i)
  !           must be set to the i-th value of the dependent variable y,
  !           for i=1,2,...,m. unchanged on exit.
  !   w     : real array of dimension at least (m). before entry, w(i)
  !           must be set to the i-th value in the set of weights. the
  !           w(i) must be strictly positive. unchanged on exit.
  !           see also further comments.
  !   xb,xe : real values. on entry xb and xe must specify the boundaries
  !           of the approximation interval. xb<=x(1), xe>=x(m).
  !           unchanged on exit.
  !   k     : integer. on entry k must specify the degree of the spline.
  !           1<=k<=5. it is recommended to use cubic splines (k=3).
  !           the user is strongly dissuaded from choosing k even,together
  !           with a small s-value. unchanged on exit.
  !   s     : real.on entry (in case iopt>=0) s must specify the smoothing
  !           factor. s >=0. unchanged on exit.
  !           for advice on the choice of s see further comments.
  !   nest  : integer. on entry nest must contain an over-estimate of the
  !           total number of knots of the spline returned, to indicate
  !           the storage space available to the routine. nest >=2*k+2.
  !           in most practical situation nest=m/2 will be sufficient.
  !           always large enough is  nest=m+k+1, the number of knots
  !           needed for interpolation (s=0). unchanged on exit.
  !   n     : integer.
  !           unless ier =10 (in case iopt >=0), n will contain the
  !           total number of knots of the spline approximation returned.
  !           if the computation mode iopt=1 is used this value of n
  !           should be left unchanged between subsequent calls.
  !           in case iopt=-1, the value of n must be specified on entry.
  !   t     : real array of dimension at least (nest).
  !           on successful exit, this array will contain the knots of the
  !           spline,i.e. the position of the interior knots t(k+2),t(k+3)
  !           ...,t(n-k-1) as well as the position of the additional knots
  !           t(1)=t(2)=...=t(k+1)=xb and t(n-k)=...=t(n)=xe needed for
  !           the b-spline representation.
  !           if the computation mode iopt=1 is used, the values of t(1),
  !           t(2),...,t(n) should be left unchanged between subsequent
  !           calls. if the computation mode iopt=-1 is used, the values
  !           t(k+2),...,t(n-k-1) must be supplied by the user, before
  !           entry. see also the restrictions (ier=10).
  !   c     : real array of dimension at least (nest).
  !           on successful exit, this array will contain the coefficients
  !           c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x)
  !   fp    : real. unless ier=10, fp contains the weighted sum of
  !           squared residuals of the spline approximation returned.
  !   wrk   : real array of dimension at least (m*(k+1)+nest*(7+3*k)).
  !           used as working space. if the computation mode iopt=1 is
  !           used, the values wrk(1),...,wrk(n) should be left unchanged
  !           between subsequent calls.
  !   lwrk  : integer. on entry,lwrk must specify the actual dimension of
  !           the array wrk as declared in the calling (sub)program.lwrk
  !           must not be too small (see wrk). unchanged on exit.
  !   iwrk  : integer array of dimension at least (nest).
  !           used as working space. if the computation mode iopt=1 is
  !           used,the values iwrk(1),...,iwrk(n) should be left unchanged
  !           between subsequent calls.
  !   ier   : integer. unless the routine detects an error, ier contains a
  !           non-positive value on exit, i.e.
  !    ier=0  : normal return. the spline returned has a residual sum of
  !             squares fp such that abs(fp-s)/s <= tol with tol a relat-
  !             ive tolerance set to 0.001 by the program.
  !    ier=-1 : normal return. the spline returned is an interpolating
  !             spline (fp=0).
  !    ier=-2 : normal return. the spline returned is the weighted least-
  !             squares polynomial of degree k. in this extreme case fp
  !             gives the upper bound fp0 for the smoothing factor s.
  !    ier=1  : error. the required storage space exceeds the available
  !             storage space, as specified by the parameter nest.
  !             probably causes : nest too small. if nest is already
  !             large (say nest > m/2), it may also indicate that s is
  !             too small
  !             the approximation returned is the weighted least-squares
  !             spline according to the knots t(1),t(2),...,t(n). (n=nest)
  !             the parameter fp gives the corresponding weighted sum of
  !             squared residuals (fp>s).
  !    ier=2  : error. a theoretically impossible result was found during
  !             the iteration process for finding a smoothing spline with
  !             fp = s. probably causes : s too small.
  !             there is an approximation returned but the corresponding
  !             weighted sum of squared residuals does not satisfy the
  !             condition abs(fp-s)/s < tol.
  !    ier=3  : error. the maximal number of iterations maxit (set to 20
  !             by the program) allowed for finding a smoothing spline
  !             with fp=s has been reached. probably causes : s too small
  !             there is an approximation returned but the corresponding
  !             weighted sum of squared residuals does not satisfy the
  !             condition abs(fp-s)/s < tol.
  !    ier=10 : error. on entry, the input data are controlled on validity
  !             the following restrictions must be satisfied.
  !             -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
  !             xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k)
  !             if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
  !                         xb<t(k+2)<t(k+3)<...<t(n-k-1)<xe
  !                       the schoenberg-whitney conditions, i.e. there
  !                       must be a subset of data points xx(j) such that
  !                         t(j) < xx(j) < t(j+k+1), j=1,2,...,n-k-1
  !             if iopt>=0: s>=0
  !                         if s=0 : nest >= m+k+1
  !             if one of these conditions is found to be violated,control
  !             is immediately repassed to the calling program. in that
  !             case there is no approximation returned.
  !
  !  further comments:
  !   by means of the parameter s, the user can control the tradeoff
  !   between closeness of fit and smoothness of fit of the approximation.
  !   if s is too large, the spline will be too smooth and signal will be
  !   lost ; if s is too small the spline will pick up too much noise. in
  !   the extreme cases the program will return an interpolating spline if
  !   s=0 and the weighted least-squares polynomial of degree k if s is
  !   very large. between these extremes, a properly chosen s will result
  !   in a good compromise between closeness of fit and smoothness of fit.
  !   to decide whether an approximation, corresponding to a certain s is
  !   satisfactory the user is highly recommended to inspect the fits
  !   graphically.
  !   recommended values for s depend on the weights w(i). if these are
  !   taken as 1/d(i) with d(i) an estimate of the standard deviation of
  !   y(i), a good s-value should be found in the range (m-sqrt(2*m),m+
  !   sqrt(2*m)). if nothing is known about the statistical error in y(i)
  !   each w(i) can be set equal to one and s determined by trial and
  !   error, taking account of the comments above. the best is then to
  !   start with a very large value of s ( to determine the least-squares
  !   polynomial and the corresponding upper bound fp0 for s) and then to
  !   progressively decrease the value of s ( say by a factor 10 in the
  !   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
  !   approximation shows more detail) to obtain closer fits.
  !   to economize the search for a good s-value the program provides with
  !   different modes of computation. at the first call of the routine, or
  !   whenever he wants to restart with the initial set of knots the user
  !   must set iopt=0.
  !   if iopt=1 the program will continue with the set of knots found at
  !   the last call of the routine. this will save a lot of computation
  !   time if curfit is called repeatedly for different values of s.
  !   the number of knots of the spline returned and their location will
  !   depend on the value of s and on the complexity of the shape of the
  !   function underlying the data. but, if the computation mode iopt=1
  !   is used, the knots returned may also depend on the s-values at
  !   previous calls (if these were smaller). therefore, if after a number
  !   of trials with different s-values and iopt=1, the user can finally
  !   accept a fit as satisfactory, it may be worthwhile for him to call
  !   curfit once more with the selected value for s but now with iopt=0.
  !   indeed, curfit may then return an approximation of the same quality
  !   of fit but with fewer knots and therefore better if data reduction
  !   is also an important objective for the user.
  !
  !  other subroutines required:
  !    fpback,fpbspl,fpchec,fpcurf,fpdisc,fpgivs,fpknot,fprati,fprota
  !
  !  references:
  !   dierckx p. : an algorithm for smoothing, differentiation and integ-
  !                ration of experimental data using spline functions,
  !                j.comp.appl.maths 1 (1975) 165-184.
  !   dierckx p. : a fast algorithm for smoothing data on a rectangular
  !                grid while using spline functions, siam j.numer.anal.
  !                19 (1982) 1286-1304.
  !   dierckx p. : an improved algorithm for curve fitting with spline
  !                functions, report tw54, dept. computer science,k.u.
  !                leuven, 1981.
  !   dierckx p. : curve and surface fitting with splines, monographs on
  !                numerical analysis, oxford university press, 1993.
  !
  !  author:
  !    p.dierckx
  !    dept. computer science, k.u. leuven
  !    celestijnenlaan 200a, b-3001 heverlee, belgium.
  !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
  !
  !  creation date : may 1979
  !  latest update : march 1987
  !
  !  ..
  !  ..scalar arguments..
  implicit none
  real(8) :: xb, xe, s, fp
  integer :: iopt, m, k, nest, n, lwrk
  integer, intent(OUT) :: ier
  !  ..array arguments..
  real(8) :: x(m), y(m), w(m), t(nest), c(nest), wrk(lwrk)
  integer :: iwrk(nest)
  !  ..local scalars..
  real(8) :: tol
  integer :: i, ia, ib, ifp, ig, iq, iz, j, k1, k2, lwest, maxit, nmin
  !  ..
  !  we set up the parameters tol and maxit
  maxit = 20
  tol = 0.001_8

  !  before starting computations a data check is made. if the input data
  !  are invalid, control is immediately repassed to the calling program.
  ier = 10
  if (k <= 0 .or. k > 5) return
  k1 = k + 1
  k2 = k1 + 1
  if (iopt < (-1) .or. iopt > 1) return
  nmin = 2 * k1
  if (m < k1 .or. nest < nmin) return
  lwest = m * k1 + nest * (7 + 3 * k)
  if (lwrk < lwest) return

  if (xb > x(1) .or. xe < x(m)) return
  do i = 2, m
    if (x(i - 1) > x(i)) return
  end do

  if (iopt >= 0) then
    if (s < 0._8) return
    if (s == 0._8 .and. nest < (m + k1)) return
  else
    if (n < nmin .or. n > nest) return
    j = n
    do i = 1, k1
      t(i) = xb
      t(j) = xe
      j = j - 1
    end do

    call fpchec(x, m, t, n, k, ier)
    if (ier /= 0) return
  end if

  ! we partition the working space and determine the spline approximation.
  ifp = 1
  iz = ifp + nest
  ia = iz + nest
  ib = ia + nest * k1
  ig = ib + nest * k2
  iq = ig + nest * k2
  call fpcurf(iopt, x, y, w, m, xb, xe, k, s, nest, tol, maxit, k1, k2, n, t, c, fp, &
    & wrk(ifp), wrk(iz), wrk(ia), wrk(ib), wrk(ig), wrk(iq), iwrk, ier)
end subroutine curfit

subroutine splder(t, n, c, k, nu, x, y, m, e, wrk, ier)
  !  subroutine splder evaluates in a number of points x(i),i=1,2,...,m
  !  the derivative of order nu of a spline s(x) of degree k,given in
  !  its b-spline representation.
  !
  !  calling sequence:
  !     call splder(t,n,c,k,nu,x,y,m,e,wrk,ier)
  !
  !  input parameters:
  !    t    : array,length n, which contains the position of the knots.
  !    n    : integer, giving the total number of knots of s(x).
  !    c    : array,length n, which contains the b-spline coefficients.
  !    k    : integer, giving the degree of s(x).
  !    nu   : integer, specifying the order of the derivative. 0<=nu<=k
  !    x    : array,length m, which contains the points where the deriv-
  !           ative of s(x) must be evaluated.
  !    m    : integer, giving the number of points where the derivative
  !           of s(x) must be evaluated
  !    e    : integer, if 0 the spline is extrapolated from the end
  !           spans for points not in the support, if 1 the spline
  !           evaluates to zero for those points, if 2 ier is set to
  !           1 and the subroutine returns, and if 3 the spline evaluates
  !           to the value of the nearest boundary point.
  !    wrk  : real array of dimension n. used as working space.
  !
  !  output parameters:
  !    y    : array,length m, giving the value of the derivative of s(x)
  !           at the different points.
  !    ier  : error flag
  !      ier = 0 : normal return
  !      ier = 1 : argument out of bounds and e == 2
  !      ier =10 : invalid input data (see restrictions)
  !
  !  restrictions:
  !    0 <= nu <= k
  !    m >= 1
  !    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
  !
  !  other subroutines required: fpbspl
  !
  !  references :
  !    de boor c : on calculating with b-splines, j. approximation theory
  !                6 (1972) 50-62.
  !    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
  !                applics 10 (1972) 134-149.
  !   dierckx p. : curve and surface fitting with splines, monographs on
  !                numerical analysis, oxford university press, 1993.
  !
  !  author :
  !    p.dierckx
  !    dept. computer science, k.u.leuven
  !    celestijnenlaan 200a, b-3001 heverlee, belgium.
  !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
  !
  !  latest update : march 1987
  !
  !++ pearu: 13 aug 2003
  !++   - disabled cliping x values to interval [min(t),max(t)]
  !++   - removed the restriction of the orderness of x values
  !++   - fixed initialization of sp to double precision value
  !
  ! JF: 2019. Converted to f20xx sintaxis and modified to use also as
  ! replacement of splev
  !  ..scalar arguments..
  implicit none
  integer, intent(IN) :: n !<
  integer, intent(IN) :: k !<
  integer, intent(IN) :: nu !<
  integer, intent(IN) :: m !<
  integer, intent(IN) :: e !<
  integer, intent(OUT) :: ier !<
  !  ..array arguments..
  real(8), dimension(n), intent(IN) :: t !<
  real(8), dimension(n), intent(IN) :: c !<
  real(8), dimension(m), intent(IN) :: x !<
  real(8), dimension(m), intent(OUT) :: y !<
  real(8), dimension(n), intent(INOUT) :: wrk !<

  !  ..local scalars..
  integer :: i, j, kk, k1, k2, l, l1, l2, nk1, nk2, nn
  real(8) :: ak, arg, fac, tb, te
  !++..
  integer :: k3
  !..++
  !  ..local arrays ..
  real(8) :: h(6)
  !  before starting computations a data check is made. if the input data
  !  are invalid control is immediately repassed to the calling program.
  ier = 10
  if (nu < 0 .or. nu > k) return
  !++..
  if (m < 1) return
  !..++
  !--  10  do i=2,m
  !--        if(x(i) < x(i-1)) return
  !--      end do
  ier = 0
  !  fetch tb and te, the boundaries of the approximation interval.
  k1 = k + 1
  k3 = k1 + 1
  nk1 = n - k1
  tb = t(k1)
  te = t(nk1 + 1)
  !  the derivative of order nu of a spline of degree k is a spline of
  !  degree k-nu,the b-spline coefficients wrk(i) of which can be found
  !  using the recurrence scheme of de boor.
  l = 1
  kk = k
  nn = n

  wrk(:nk1) = c(:nk1)

  if (nu > 0) then
    nk2 = nk1
    do j = 1, nu
      ak = kk
      nk2 = nk2 - 1
      l1 = l
      do i = 1, nk2
        l1 = l1 + 1
        l2 = l1 + kk
        fac = t(l2) - t(l1)
        IF (fac > 0._8) wrk(i) = ak * (wrk(i + 1) - wrk(i)) / fac
      end do

      l = l + 1
      kk = kk - 1
    end do
  end if

  if (kk == 0) then !  if nu=k the derivative is a piecewise constant function
    j = 1
  else
    l = k1
    k2 = k1 - nu
  end if

  !  main loop for the different points.
  do i = 1, m
    !  fetch a new x-value arg.
    arg = x(i)
    !  check if arg is in the support
    if (arg < tb .or. arg > te) then
      if (e == 3) then
        if (arg < tb) then
          arg = tb
        else
          arg = te
        endif
      else if (e == 2) then
        ier = 1
        return
      else if (e == 1) then
        y(i) = 0
        cycle
      endif
    endif
    !  search for knot interval t(l) <= arg < t(l+1)
    !JF: This should be quite efficient if x(i) are ordered.
    do while (arg < t(l) .and. l /= k1)
      l = l - 1
      j = j - 1
    end do
    do while (arg >= t(l + 1) .and. l /= nk1)
      l = l + 1
      j = j + 1
    end do
    !  evaluate the non-zero b-splines of degree k-nu at arg.
    if (kk == 0) then        !  nu=k: the derivative is a piecewise function
      y(i) = wrk(j)
    else
      call fpbspl(t, n, kk, arg, l, h)
      y(i) = sum(wrk(l - k1 + 1:l - nu) * h(:k2))
    end if
  end do
end subroutine splder

subroutine sproot(t, n, c, zeros, mest, m, ier)
  !  subroutine sproot finds the zeros of a cubic spline s(x),which is
  !  given in its normalized b-spline representation.
  !
  !  calling sequence:
  !     call sproot(t,n,c,zeros,mest,m,ier)
  !
  !  input parameters:
  !    t    : real array,length n, containing the knots of s(x).
  !    n    : integer, containing the number of knots.  n>=8
  !    c    : real array,length n, containing the b-spline coefficients.
  !    mest : integer, specifying the dimension of array zero.
  !
  !  output parameters:
  !    zeros : real array,length mest, containing the zeros of s(x).
  !    m    : integer,giving the number of zeros.
  !    ier  : error flag:
  !      ier = 0: normal return.
  !      ier = 1: the number of zeros exceeds mest.
  !      ier =10: invalid input data (see restrictions).
  !
  !  other subroutines required: fpcuro
  !
  !  restrictions:
  !    1) n>= 8.
  !    2) t(4) < t(5) < ... < t(n-4) < t(n-3).
  !       t(1) <= t(2) <= t(3) <= t(4)
  !       t(n-3) <= t(n-2) <= t(n-1) <= t(n)
  !
  !  author :
  !    p.dierckx
  !    dept. computer science, k.u.leuven
  !    celestijnenlaan 200a, b-3001 heverlee, belgium.
  !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
  !
  !  latest update : march 1987
  !
  ! ..
  ! ..scalar arguments..
  implicit none
  integer, intent(IN) :: n !< n = size(t)
  real(8), dimension(n), intent(IN) :: t !< Knots
  real(8), dimension(n), intent(IN) :: c !< Coefficients
  integer, intent(IN) :: mest !<
  real(8), dimension(mest), intent(OUT) :: zeros !<
  integer, intent(OUT) :: m !<
  integer, intent(OUT) :: ier !<
  !  ..local scalars..
  integer :: i, j, j1, l, n4
  real(8) :: ah, a0, a1, a2, a3, bh, b0, b1, c1, c2, c3, c4, c5, d4, d5, h1, h2,&
    &three, two, t1, t2, t3, t4, t5, zz
  logical z0, z1, z2, z3, z4, nz0, nz1, nz2, nz3, nz4
  !  ..local array..
  real(8) :: y(3)
  !  ..
  !  set some constants
  two = 2._8
  three = 3._8
  !  before starting computations a data check is made. if the input data
  !  are invalid, control is immediately repassed to the calling program.
  n4 = n - 4
  ier = 10
  if (n < 8) return
  j = n
  do i = 1, 3
    if (t(i) > t(i + 1)) return
    if (t(j) < t(j - 1)) return
    j = j - 1
  end do

  do i = 4, n4
    if (t(i) >= t(i + 1)) return
  end do

  !  the problem considered reduces to finding the zeros of the cubic
  !  polynomials pl(x) which define the cubic spline in each knot
  !  interval t(l)<=x<=t(l+1). a zero of pl(x) is also a zero of s(x) on
  !  the condition that it belongs to the knot interval.
  !  the cubic polynomial pl(x) is determined by computing s(t(l)),
  !  s'(t(l)),s(t(l+1)) and s'(t(l+1)). in fact we only have to compute
  !  s(t(l+1)) and s'(t(l+1)); because of the continuity conditions of
  !  splines and their derivatives, the value of s(t(l)) and s'(t(l))
  !  is already known from the foregoing knot interval.
  ier = 0
  !  evaluate some constants for the first knot interval
  h1 = t(4) - t(3)
  h2 = t(5) - t(4)
  t1 = t(4) - t(2)
  t2 = t(5) - t(3)
  t3 = t(6) - t(4)
  t4 = t(5) - t(2)
  t5 = t(6) - t(3)
  !  calculate a0 = s(t(4)) and ah = s'(t(4)).
  c1 = c(1)
  c2 = c(2)
  c3 = c(3)
  c4 = (c2 - c1) / t4
  c5 = (c3 - c2) / t5
  d4 = (h2 * c1 + t1 * c2) / t4
  d5 = (t3 * c2 + h1 * c3) / t5
  a0 = (h2 * d4 + h1 * d5) / t2
  ah = three * (h2 * c4 + h1 * c5) / t2
  z1 = .true.
  if (ah < 0._8) z1 = .false.
  nz1 = .not. z1
  m = 0
  !  main loop for the different knot intervals.
  do l = 4, n4
    !  evaluate some constants for the knot interval t(l) <= x <= t(l+1).
    h1 = h2
    h2 = t(l + 2) - t(l + 1)
    t1 = t2
    t2 = t3
    t3 = t(l + 3) - t(l + 1)
    t4 = t5
    t5 = t(l + 3) - t(l)
    !  find a0 = s(t(l)), ah = s'(t(l)), b0 = s(t(l+1)) and bh = s'(t(l+1)).
    c1 = c2
    c2 = c3
    c3 = c(l)
    c4 = c5
    c5 = (c3 - c2) / t5
    d4 = (h2 * c1 + t1 * c2) / t4
    d5 = (h1 * c3 + t3 * c2) / t5
    b0 = (h2 * d4 + h1 * d5) / t2
    bh = three * (h2 * c4 + h1 * c5) / t2
    !  calculate the coefficients a0,a1,a2 and a3 of the cubic polynomial
    !  pl(x) = ql(y) = a0+a1*y+a2*y**2+a3*y**3 ; y = (x-t(l))/(t(l+1)-t(l)).
    a1 = ah * h1
    b1 = bh * h1
    a2 = three * (b0 - a0) - b1 - two * a1
    a3 = two * (a0 - b0) + b1 + a1
    !  test whether or not pl(x) could have a zero in the range
    !  t(l) <= x <= t(l+1).
    z3 = .true.
    if (b1 < 0._8) z3 = .false.
    nz3 = .not. z3
    if (a0 * b0 <= 0._8) go to 100
    z0 = .true.
    if (a0 < 0._8) z0 = .false.
    nz0 = .not. z0
    z2 = .true.
    if (a2 < 0.) z2 = .false.
    nz2 = .not. z2
    z4 = .true.
    if (3.0d0 * a3 + a2 < 0._8) z4 = .false.
    nz4 = .not. z4
    if (.not. ((z0 .and. (nz1 .and. (z3 .or. z2 .and. nz4) .or. nz2 .and.&
      &z3 .and. z4) .or. nz0 .and. (z1 .and. (nz3 .or. nz2 .and. z4) .or. z2 .and.&
      &nz3 .and. nz4)))) go to 200
    !  find the zeros of ql(y).
100 call fpcuro(a3, a2, a1, a0, y, j)
    if (j == 0) go to 200
    !  find which zeros of pl(x) are zeros of s(x).
    do i = 1, j
      if (y(i) < 0.0_8 .or. y(i) > 1.0_8) go to 150
      !  test whether the number of zeros of s(x) exceeds mest.
      if (m >= mest) then
        ier = 1
        return
      end if

      m = m + 1
      zeros(m) = t(l) + h1 * y(i)
150   continue
    end do

200 a0 = b0
    ah = bh
    z1 = z3
    nz1 = nz3
  end do

  !  the zeros of s(x) are arranged in increasing order.
  if (m < 2) return
  do i = 2, m
    j = i
350 j1 = j - 1
    if (j1 == 0) go to 400
    if (zeros(j) >= zeros(j1)) go to 400
    zz = zeros(j)
    zeros(j) = zeros(j1)
    zeros(j1) = zz
    j = j1
    go to 350
400 continue
  end do

  j = m
  m = 1
  do i = 2, j
    if (zeros(i) /= zeros(m)) then
      m = m + 1
      zeros(m) = zeros(i)
    end if
  end do

end subroutine sproot

! !!!!!!!!!!!!!!!!!!!!  Auxiliary routines !!!!!!!!!!!!!!!!!!!!!

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fpback(a, z, n, k, c, nest)
  !  subroutine fpback calculates the solution of the system of
  !  equations a*c = z with a a n x n upper triangular matrix
  !  of bandwidth k.
  !  ..
  !  ..scalar arguments..
  integer :: n, k, nest
  !  ..array arguments..
  real(8) :: a(nest, k), z(n), c(n)
  !  ..local scalars..
  real(8) :: store
  integer :: i, i1, j, k1, l, m
  !  ..

  k1 = k - 1
  c(n) = z(n) / a(n, 1)
  i = n - 1
  if (i == 0) return
  do j = 2, n
    store = z(i)
    i1 = k1
    if (j <= k1) i1 = j - 1
    m = i
    do l = 1, i1
      m = m + 1
      store = store - c(m) * a(i, l + 1)
    end do
    c(i) = store / a(i, 1)
    i = i - 1
  end do
end subroutine fpback
subroutine fpbacp(a, b, z, n, k, c, k1, nest)
  !  subroutine fpbacp calculates the solution of the system of equations
  !  g * c = z  with g  a n x n upper triangular matrix of the form
  !            ! a '   !
  !        g = !   ' b !
  !            ! 0 '   !
  !  with b a n x k matrix and a a (n-k) x (n-k) upper triangular
  !  matrix of bandwidth k1.
  !  ..
  !  ..scalar arguments..
  integer :: n, k, k1, nest
  !  ..array arguments..
  real(8) :: a(nest, k1), b(nest, k), z(n), c(n)
  !  ..local scalars..
  integer :: i, i1, j, l, l0, l1, n2
  real(8) :: store
  !  ..
  n2 = n - k
  l = n
  do i = 1, k
    store = z(l)
    j = k + 2 - i
    if (i > 1) then
      l0 = l
      do l1 = j, k
        l0 = l0 + 1
        store = store - c(l0) * b(l, l1)
      end do
    end if

    c(l) = store / b(l, j - 1)
    l = l - 1
    if (l == 0) return
  end do

  do i = 1, n2
    store = z(i)
    l = n2
    do j = 1, k
      l = l + 1
      store = store - c(l) * b(i, j)
    end do

    c(i) = store
  end do

  i = n2
  c(i) = c(i) / a(i, 1)

  if (i == 1) return

  do j = 2, n2
    i = i - 1
    store = c(i)
    i1 = k
    if (j <= k) i1 = j - 1
    l = i
    do l0 = 1, i1
      l = l + 1
      store = store - c(l) * a(i, l0 + 1)
    end do
    c(i) = store / a(i, 1)
  end do

end subroutine fpbacp
subroutine fpbspl(t, n, k, x, l, h)
  !  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
  !  degree k at t(l) <= x < t(l+1) using the stable recurrence
  !  relation of de boor and cox.
  !  Travis Oliphant  2007
  !    changed so that weighting of 0 is used when knots with
  !      multiplicity are present.
  !    Also, notice that l+k <= n and 1 <= l+1-k
  !      or else the routine will be accessing memory outside t
  !      Thus it is imperative that that k <= l <= n-k but this
  !      is not checked.
  !  ..
  !  ..scalar arguments..
  real(8), intent(IN) :: x
  integer, intent(IN) :: n, k, l
  !  ..array arguments..
  real(8), intent(IN), dimension(n):: t
  real(8), intent(OUT), dimension(:) :: h
  !  ..local scalars..
  real(8) :: f
  integer :: i, j, li, lj
  !  ..local arrays..
  real(8) :: hh(19)
  !  ..

  h(1) = 1._8
  do j = 1, k
    hh(:j) = h(:j)
    h(1) = 0.0_8
    do i = 1, j
      li = l + i
      lj = li - j
      if (t(li) == t(lj)) then
        h(i + 1) = 0.0_8
      else
        f = hh(i) / (t(li) - t(lj))
        h(i) = h(i) + f * (t(li) - x)
        h(i + 1) = f * (x - t(lj))
      end if
    end do
  end do

end subroutine fpbspl
subroutine fpchec(x, m, t, n, k, ier)
  !  subroutine fpchec verifies the number and the position of the knots
  !  t(j),j=1,2,...,n of a spline of degree k, in relation to the number
  !  and the position of the data points x(i),i=1,2,...,m. if all of the
  !  following conditions are fulfilled, the error parameter ier is set
  !  to zero. if one of the conditions is violated ier is set to ten.
  !      1) k+1 <= n-k-1 <= m
  !      2) t(1) <= t(2) <= ... <= t(k+1)
  !         t(n-k) <= t(n-k+1) <= ... <= t(n)
  !      3) t(k+1) < t(k+2) < ... < t(n-k)
  !      4) t(k+1) <= x(i) <= t(n-k)
  !      5) the conditions specified by schoenberg and whitney must hold
  !         for at least one subset of data points, i.e. there must be a
  !         subset of data points y(j) such that
  !             t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
  !  ..
  !  ..scalar arguments..
  integer :: m, n, k, ier
  !  ..array arguments..
  real(8) :: x(m), t(n)
  !  ..local scalars..
  integer :: i, j, k1, k2, l, nk1, nk2, nk3
  real(8) :: tj, tl
  !  ..
  k1 = k + 1
  k2 = k1 + 1
  nk1 = n - k1
  nk2 = nk1 + 1
  ier = 10
  !  check condition no 1
  if (nk1 < k1 .or. nk1 > m) return
  !  check condition no 2
  j = n
  do i = 1, k
    if (t(i) > t(i + 1)) return
    if (t(j) < t(j - 1)) return
    j = j - 1
  end do
  !  check condition no 3
  do i = k2, nk2
    if (t(i) <= t(i - 1)) return
  end do
  !  check condition no 4
  if (x(1) < t(k1) .or. x(m) > t(nk2)) return
  !  check condition no 5
  if (x(1) >= t(k2) .or. x(m) <= t(nk1)) return
  i = 1
  l = k2
  nk3 = nk1 - 1
  if (nk3 < 2) then
    ier = 0
    return
  end if

  do j = 2, nk3
    tj = t(j)
    l = l + 1
    tl = t(l)
40  i = i + 1
    if (i >= m) return
    if (x(i) <= tj) go to 40
    if (x(i) >= tl) return
  end do

end subroutine fpchec
subroutine fpchep(x, m, t, n, k, ier)
  !  subroutine fpchep verifies the number and the position of the knots
  !  t(j),j=1,2,...,n of a periodic spline of degree k, in relation to
  !  the number and the position of the data points x(i),i=1,2,...,m.
  !  if all of the following conditions are fulfilled, ier is set
  !  to zero. if one of the conditions is violated ier is set to ten.
  !      1) k+1 <= n-k-1 <= m+k-1
  !      2) t(1) <= t(2) <= ... <= t(k+1)
  !         t(n-k) <= t(n-k+1) <= ... <= t(n)
  !      3) t(k+1) < t(k+2) < ... < t(n-k)
  !      4) t(k+1) <= x(i) <= t(n-k)
  !      5) the conditions specified by schoenberg and whitney must hold
  !         for at least one subset of data points, i.e. there must be a
  !         subset of data points y(j) such that
  !             t(j) < y(j) < t(j+k+1), j=k+1,...,n-k-1
  !  ..
  !  ..scalar arguments..
  integer :: m, n, k, ier
  !  ..array arguments..
  real(8) :: x(m), t(n)
  !  ..local scalars..
  integer :: i, i1, i2, j, j1, k1, k2, l, l1, l2, mm, m1, nk1, nk2
  real(8) :: per, tj, tl, xi
  !  ..
  k1 = k + 1
  k2 = k1 + 1
  nk1 = n - k1
  nk2 = nk1 + 1
  m1 = m - 1
  ier = 10
  !  check condition no 1
  if (nk1 < k1 .or. n > m + 2 * k) return
  !  check condition no 2
  j = n
  do i = 1, k
    if (t(i) > t(i + 1)) return
    if (t(j) < t(j - 1)) return
    j = j - 1
  end do

  !  check condition no 3
  do i = k2, nk2
    if (t(i) <= t(i - 1)) return
  end do

  !  check condition no 4
  if (x(1) < t(k1) .or. x(m) > t(nk2)) return
  !  check condition no 5
  l1 = k1
  l2 = 1
  do l = 1, m
    xi = x(l)
40  if (xi < t(l1 + 1) .or. l == nk1) cycle
    l1 = l1 + 1
    l2 = l2 + 1
    if (l2 > k1) go to 60
    go to 40
  end do

  l = m
60 per = t(nk2) - t(k1)
  do i1 = 2, l
    i = i1 - 1
    mm = i + m1
    do j = k1, nk1
      tj = t(j)
      j1 = j + k1
      tl = t(j1)
70    i = i + 1
      if (i > mm) exit
      i2 = i - m1
      if (i2 <= 0) go to 80
      go to 90
80    xi = x(i)
      go to 100
90    xi = x(i2) + per
100   if (xi <= tj) go to 70
      if (xi >= tl) exit
    end do

    ier = 0
    return
  end do

end subroutine fpchep
subroutine fpclos(iopt, idim, m, u, mx, x, w, k, s, nest, tol, maxit, k1, k2,&
  &n, t, nc, c, fp, fpint, z, a1, a2, b, g1, g2, q, nrdata, ier)
  !  ..
  !  ..scalar arguments..

  real(8) :: s, tol, fp
  integer :: iopt, idim, m, mx, k, nest, maxit, k1, k2, n, nc, ier
  !  ..array arguments..
  real(8) :: u(m), x(mx), w(m), t(nest), c(nc), fpint(nest), z(nc), a1(nest, k1),&
    & a2(nest, k), b(nest, k2), g1(nest, k2), g2(nest, k1), q(m, k1)
  integer :: nrdata(nest)
  !  ..local scalars..
  real(8) :: acc, cos, d1, fac, fpart, fpms, fpold, fp0, f1, f2, f3, p, per, pinv, piv,&
    &p1, p2, p3, sin, store, term, ui, wi, rn, one, con1, con4, con9, half
  integer :: i, ich1, ich3, ij, ik, it, iter, i1, i2, i3, j, jj, jk, jper, j1, j2, kk,&
    &kk1, k3, l, l0, l1, l5, mm, m1, new, nk1, nk2, nmax, nmin, nplus, npl1,&
    &nrint, n10, n11, n7, n8
  !  ..local arrays..
  real(8) :: h(6), h1(7), h2(6), xi(10)
  !  ..function references..
  ! real(8) :: abs, fprati
  ! integer :: max0, min0
  !  ..subroutine references..
  !    fpbacp,fpbspl,fpgivs,fpdisc,fpknot,fprota
  !  ..
  !  set constants
  one = 0.1e+01_8
  con1 = 0.1e0_8
  con9 = 0.9e0_8
  con4 = 0.4e-01_8
  half = 0.5e0_8
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  part 1: determination of the number of knots and their position     c
  !  **************************************************************      c
  !  given a set of knots we compute the least-squares closed curve      c
  !  sinf(u). if the sum f(p=inf) <= s we accept the choice of knots.    c
  !  if iopt=-1 sinf(u) is the requested curve                           c
  !  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
  !    if fp <=s we will continue with the current set of knots.         c
  !    if fp > s we will increase the number of knots and compute the    c
  !       corresponding least-squares curve until finally fp<=s.         c
  !  the initial choice of knots depends on the value of s and iopt.     c
  !    if s=0 we have spline interpolation; in that case the number of   c
  !    knots equals nmax = m+2*k.                                        c
  !    if s > 0 and                                                      c
  !      iopt=0 we first compute the least-squares polynomial curve of   c
  !      degree k; n = nmin = 2*k+2. since s(u) must be periodic we      c
  !      find that s(u) reduces to a fixed point.                        c
  !      iopt=1 we start with the set of knots found at the last         c
  !      call of the routine, except for the case that s > fp0; then     c
  !      we compute directly the least-squares polynomial curve.         c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  m1 = m - 1
  kk = k
  kk1 = k1
  k3 = 3 * k + 1
  nmin = 2 * k1
  !  determine the length of the period of the splines.
  per = u(m) - u(1)
  if (iopt < 0) go to 50
  !  calculation of acc, the absolute tolerance for the root of f(p)=s.
  acc = tol * s
  !  determine nmax, the number of knots for periodic spline interpolation
  nmax = m + 2 * k
  if (s > 0. .or. nmax == nmin) go to 30
  !  if s=0, s(u) is an interpolating curve.
  n = nmax
  !  test whether the required storage space exceeds the available one.
  if (n > nest) then
    ier = 1
    return
  end if

  !  find the position of the interior knots in case of interpolation.
5 if ((k / 2) * 2 == k) go to 20
  do i = 2, m1
    j = i + k
    t(j) = u(i)
  end do

  if (s > 0.) go to 50
  kk = k - 1
  kk1 = k
  if (kk > 0) go to 50
  t(1) = t(m) - per
  t(2) = u(1)
  t(m + 1) = u(m)
  t(m + 2) = t(3) + per
  jj = 0
  do i = 1, m1
    j = i
    do j1 = 1, idim
      jj = jj + 1
      c(j) = x(jj)
      j = j + n
    end do
  end do

  jj = 1
  j = m
  do j1 = 1, idim
    c(j) = c(jj)
    j = j + n
    jj = jj + n
  end do

  fp = 0.
  fpint(n) = fp0
  fpint(n - 1) = 0.
  nrdata(n) = 0
  ier = -1
  return

20 do i = 2, m1
    j = i + k
    t(j) = (u(i) + u(i - 1)) * half
  end do

  go to 50
  !  if s > 0 our initial choice depends on the value of iopt.
  !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
  !  polynomial curve. (i.e. a constant point).
  !  if iopt=1 and fp0>s we start computing the least-squares closed
  !  curve according the set of knots found at the last call of the
  !  routine.
30 if (iopt == 0) go to 35
  if (n == nmin) go to 35
  fp0 = fpint(n)
  fpold = fpint(n - 1)
  nplus = nrdata(n)
  if (fp0 > s) go to 50
  !  the case that s(u) is a fixed point is treated separetely.
  !  fp0 denotes the corresponding sum of squared residuals.
35 fp0 = 0.
  d1 = 0.
  z(:idim) = 0._8
  jj = 0
  do it = 1, m1
    wi = w(it)
    call fpgivs(wi, d1, cos, sin)
    do j = 1, idim
      jj = jj + 1
      fac = wi * x(jj)
      call fprota(cos, sin, fac, z(j))
      fp0 = fp0 + fac**2
    end do
  end do

  z(:idim) = z(:idim) / d1

  !  test whether that fixed point is a solution of our problem.
  fpms = fp0 - s
  if (fpms < acc .or. nmax == nmin) go to 640
  fpold = fp0
  !  test whether the required storage space exceeds the available one.
  if (n >= nest) then
    ier = 1
    return
  end if

  !  start computing the least-squares closed curve with one
  !  interior knot.
  nplus = 1
  n = nmin + 1
  mm = (m + 1) / 2
  t(k2) = u(mm)
  nrdata(1) = mm - 2
  nrdata(2) = m1 - mm
  !  main loop for the different sets of knots. m is a save upper
  !  bound for the number of trials.
50 do iter = 1, m
    !  find nrint, the number of knot intervals.
    nrint = n - nmin + 1
    !  find the position of the additional knots which are needed for
    !  the b-spline representation of s(u). if we take
    !      t(k+1) = u(1), t(n-k) = u(m)
    !      t(k+1-j) = t(n-k-j) - per, j=1,2,...k
    !      t(n-k+j) = t(k+1+j) + per, j=1,2,...k
    !  then s(u) will be a smooth closed curve if the b-spline
    !  coefficients satisfy the following conditions
    !      c((i-1)*n+n7+j) = c((i-1)*n+j), j=1,...k,i=1,2,...,idim (**)
    !  with n7=n-2*k-1.
    t(k1) = u(1)
    nk1 = n - k1
    nk2 = nk1 + 1
    t(nk2) = u(m)
    do j = 1, k
      i1 = nk2 + j
      i2 = nk2 - j
      j1 = k1 + j
      j2 = k1 - j
      t(i1) = t(j1) + per
      t(j2) = t(i2) - per
    end do

    !  compute the b-spline coefficients of the least-squares closed curve
    !  sinf(u). the observation matrix a is built up row by row while
    !  taking into account condition (**) and is reduced to triangular
    !  form by givens transformations .
    !  at the same time fp=f(p=inf) is computed.
    !  the n7 x n7 triangularised upper matrix a has the form
    !            ! a1 '    !
    !        a = !    ' a2 !
    !            ! 0  '    !
    !  with a2 a n7 x k matrix and a1 a n10 x n10 upper triangular
    !  matrix of bandwidth k+1 ( n10 = n7-k).
    !  initialization.

    z(:nc) = 0._8
    a1(:nk1, :kk1) = 0._8

    n7 = nk1 - k
    n10 = n7 - kk
    jper = 0
    fp = 0.
    l = k1
    jj = 0
    item: do it = 1, m1
      !  fetch the current data point u(it),x(it)
      ui = u(it)
      wi = w(it)
      do j = 1, idim
        jj = jj + 1
        xi(j) = x(jj) * wi
      end do

      !  search for knot interval t(l) <= ui < t(l+1).
80    if (ui < t(l + 1)) go to 85
      l = l + 1
      go to 80
      !  evaluate the (k+1) non-zero b-splines at ui and store them in q.
85    call fpbspl(t, n, k, ui, l, h)
      do i = 1, k1
        q(it, i) = h(i)
        h(i) = h(i) * wi
      end do

      l5 = l - k1
      !  test whether the b-splines nj,k+1(u),j=1+n7,...nk1 are all zero at ui
      if (l5 < n10) go to 285
      if (jper /= 0) go to 160
      !  initialize the matrix a2.
      a2(:n7, :kk) = 0._8
      jk = n10 + 1
      do i = 1, kk
        ik = jk
        do j = 1, kk1
          if (ik <= 0) cycle
          a2(ik, i) = a1(ik, j)
          ik = ik - 1
        end do

        jk = jk + 1
      end do

      jper = 1
      !  if one of the b-splines nj,k+1(u),j=n7+1,...nk1 is not zero at ui
      !  we take account of condition (**) for setting up the new row
      !  of the observation matrix a. this row is stored in the arrays h1
      !  (the part with respect to a1) and h2 (the part with
      !  respect to a2).
160   h1(:kk) = 0._8
      h2(:kk) = 0._8
      h1(kk1) = 0.
      j = l5 - n10
      do i = 1, kk1
        j = j + 1
        l0 = j
        do
180       l1 = l0 - kk
          if (l1 <= 0) go to 200
          if (l1 <= n10) go to 190
          l0 = l1 - n10
          go to 180
        end do

190     h1(l1) = h(i)
        cycle
200     h2(l0) = h2(l0) + h(i)

      end do

      !  rotate the new row of the observation matrix into triangle
      !  by givens transformations.
      if (n10 <= 0) go to 250
      !  rotation with the rows 1,2,...n10 of matrix a.
      do j = 1, n10
        piv = h1(1)
        if (piv /= 0.) go to 214

        do i = 1, kk
          h1(i) = h1(i + 1)
        end do
        h1(kk1) = 0.
        cycle ! go to 240
        !  calculate the parameters of the givens transformation.
214     call fpgivs(piv, a1(j, 1), cos, sin)
        !  transformation to the right hand side.
        j1 = j
        do j2 = 1, idim
          call fprota(cos, sin, xi(j2), z(j1))
          j1 = j1 + n
        end do

        !  transformations to the left hand side with respect to a2.
        do i = 1, kk
          call fprota(cos, sin, h2(i), a2(j, i))
        end do

        if (j == n10) go to 250
        i2 = min0(n10 - j, kk)
        !  transformations to the left hand side with respect to a1.
        do i = 1, i2
          i1 = i + 1
          call fprota(cos, sin, h1(i1), a1(j, i1))
          h1(i) = h1(i1)
        end do

        h1(i1) = 0.

      end do

      !  rotation with the rows n10+1,...n7 of matrix a.
250   do j = 1, kk
        ij = n10 + j
        if (ij <= 0) cycle
        piv = h2(j)
        if (piv == 0.) cycle
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, a2(ij, j), cos, sin)
        !  transformations to right hand side.
        j1 = ij
        do j2 = 1, idim
          call fprota(cos, sin, xi(j2), z(j1))
          j1 = j1 + n
        end do

        if (j == kk) go to 280
        j1 = j + 1
        !  transformations to left hand side.
        do i = j1, kk
          call fprota(cos, sin, h2(i), a2(ij, i))
        end do

      end do

      !  add contribution of this row to the sum of squares of residual
      !  right hand sides.
280   do j2 = 1, idim
        fp = fp + xi(j2)**2
      end do

      cycle item
      !  rotation of the new row of the observation matrix into
      !  triangle in case the b-splines nj,k+1(u),j=n7+1,...n-k-1 are all zero
      !  at ui.
285   j = l5
      do i = 1, kk1
        j = j + 1
        piv = h(i)
        if (piv == 0.) cycle
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, a1(j, 1), cos, sin)
        !  transformations to right hand side.
        j1 = j
        do j2 = 1, idim
          call fprota(cos, sin, xi(j2), z(j1))
          j1 = j1 + n
        end do

        if (i == kk1) exit
        i2 = 1
        i3 = i + 1
        !  transformations to left hand side.
        do i1 = i3, kk1
          i2 = i2 + 1
          call fprota(cos, sin, h(i1), a1(j, i2))
        end do

      end do

      !  add contribution of this row to the sum of squares of residual
      !  right hand sides.
      do j2 = 1, idim
        fp = fp + xi(j2)**2
      end do

    end do item

    fpint(n) = fp0
    fpint(n - 1) = fpold
    nrdata(n) = nplus
    !  backward substitution to obtain the b-spline coefficients .
    j1 = 1
    do j2 = 1, idim
      call fpbacp(a1, a2, z(j1), n7, kk, c(j1), kk1, nest)
      j1 = j1 + n
    end do

    !  calculate from condition (**) the remaining coefficients.
    do i = 1, k
      j1 = i
      do j = 1, idim
        j2 = j1 + n7
        c(j2) = c(j1)
        j1 = j1 + n
      end do
    end do

    if (iopt < 0) return
    !  test whether the approximation sinf(u) is an acceptable solution.
    fpms = fp - s
    if (abs(fpms) < acc) return
    !  if f(p=inf) < s accept the choice of knots.
    if (fpms < 0.) go to 350
    !  if n=nmax, sinf(u) is an interpolating curve.
    if (n == nmax) then
      ier = -1
      return
    end if

    !  increase the number of knots.
    !  if n=nest we cannot increase the number of knots because of the
    !  storage capacity limitation.
    if (n == nest) then
      ier = 1
      return
    end if

    !  determine the number of knots nplus we are going to add.
    npl1 = nplus * 2
    rn = nplus
    if (fpold - fp > acc) npl1 = int(rn * fpms / (fpold - fp))
    nplus = min(nplus * 2, max(npl1, nplus / 2, 1))
    fpold = fp
    !  compute the sum of squared residuals for each knot interval
    !  t(j+k) <= ui <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
    fpart = 0.
    i = 1
    l = k1
    jj = 0
    do it = 1, m1
      if (u(it) < t(l)) go to 300
      new = 1
      l = l + 1
300   term = 0.
      l0 = l - k2
      do j2 = 1, idim
        fac = 0.
        j1 = l0
        do j = 1, k1
          j1 = j1 + 1
          fac = fac + c(j1) * q(it, j)
        end do

        jj = jj + 1
        term = term + (w(it) * (fac - x(jj)))**2
        l0 = l0 + n
      end do

      fpart = fpart + term
      if (new == 0) cycle
      if (l > k2) go to 315
      fpint(nrint) = term
      new = 0
      cycle
315   store = term * half
      fpint(i) = fpart - store
      i = i + 1
      fpart = store
      new = 0
    end do

    fpint(nrint) = fpint(nrint) + fpart
    knots: do l = 1, nplus
      !  add a new knot
      call fpknot(u, m, t, n, fpint, nrdata, nrint, nest, 1)
      !  if n=nmax we locate the knots as for interpolation
      if (n == nmax) go to 5
      !  test whether we cannot further increase the number of knots.
      if (n == nest) exit knots
    end do knots

    !  restart the computations with the new set of knots.
  end do

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  part 2: determination of the smoothing closed curve sp(u).          c
  !  **********************************************************          c
  !  we have determined the number of knots and their position.          c
  !  we now compute the b-spline coefficients of the smoothing curve     c
  !  sp(u). the observation matrix a is extended by the rows of matrix   c
  !  b expressing that the kth derivative discontinuities of sp(u) at    c
  !  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
  !  ponding weights of these additional rows are set to 1/p.            c
  !  iteratively we then have to determine the value of p such that f(p),c
  !  the sum of squared residuals be = s. we already know that the least-c
  !  squares polynomial curve corresponds to p=0, and that the least-    c
  !  squares periodic spline curve corresponds to p=infinity. the        c
  !  iteration process which is proposed here, makes use of rational     c
  !  interpolation. since f(p) is a convex and strictly decreasing       c
  !  function of p, it can be approximated by a rational function        c
  !  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
  !  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
  !  to calculate the new value of p such that r(p)=s. convergence is    c
  !  guaranteed by taking f1>0 and f3<0.                                 c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  evaluate the discontinuity jump of the kth derivative of the
  !  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
350 call fpdisc(t, n, k2, b, nest)
  !  initial value for p.
  p1 = 0.
  f1 = fp0 - s
  p3 = -one
  f3 = fpms
  n11 = n10 - 1
  n8 = n7 - 1
  p = 0.
  l = n7
  do i = 1, k
    j = k + 1 - i
    p = p + a2(l, j)
    l = l - 1
    if (l == 0) go to 356
  end do

  p = p + sum(a1(:n10, 1))
  ! do i = 1, n10
  !   p = p + a1(i, 1)
  ! end do

356 rn = n7
  p = rn / p
  ich1 = 0
  ich3 = 0
  !  iteration process to find the root of f(p) = s.
  do iter = 1, maxit
    !  form the matrix g  as the matrix a extended by the rows of matrix b.
    !  the rows of matrix b with weight 1/p are rotated into
    !  the triangularised observation matrix a.
    !  after triangularisation our n7 x n7 matrix g takes the form
    !            ! g1 '    !
    !        g = !    ' g2 !
    !            ! 0  '    !
    !  with g2 a n7 x (k+1) matrix and g1 a n11 x n11 upper triangular
    !  matrix of bandwidth k+2. ( n11 = n7-k-1)
    pinv = one / p
    !  store matrix a into g
    c(:nc) = z(:nc)
    do i = 1, n7
      g1(i, k1) = a1(i, k1)
      g1(i, k2) = 0.
      g2(i, 1) = 0.
      do j = 1, k
        g1(i, j) = a1(i, j)
        g2(i, j + 1) = a2(i, j)
      end do
    end do

    l = n10
    do j = 1, k1
      if (l <= 0) exit
      g2(l, 1) = a1(l, j)
      l = l - 1
    end do

    do it = 1, n8
      !  fetch a new row of matrix b and store it in the arrays h1 (the part
      !  with respect to g1) and h2 (the part with respect to g2).
      xi(:idim) = 0._8
      h1(:k1) = 0._8
      h2(:k1) = 0._8
      h1(k2) = 0.
      if (it > n11) go to 420
      l = it
      l0 = it
      do j = 1, k2
        if (l0 == n10) go to 400
        h1(j) = b(it, j) * pinv
        l0 = l0 + 1
      end do
      go to 470

400   l0 = 1
      do l1 = j, k2
        h2(l0) = b(it, l1) * pinv
        l0 = l0 + 1
      end do
      go to 470

420   l = 1
      i = it - n10
      do j = 1, k2
        i = i + 1
        l0 = i
430     l1 = l0 - k1
        if (l1 <= 0) go to 450
        if (l1 <= n11) go to 440
        l0 = l1 - n11
        go to 430
440     h1(l1) = b(it, j) * pinv
        cycle
450     h2(l0) = h2(l0) + b(it, j) * pinv
      end do

      if (n11 <= 0) go to 510
      !  rotate this row into triangle by givens transformations
      !  rotation with the rows l,l+1,...n11.
470   do j = l, n11
        piv = h1(1)
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, g1(j, 1), cos, sin)
        !  transformation to right hand side.
        j1 = j
        do j2 = 1, idim
          call fprota(cos, sin, xi(j2), c(j1))
          j1 = j1 + n
        end do

        !  transformation to the left hand side with respect to g2.
        do i = 1, k1
          call fprota(cos, sin, h2(i), g2(j, i))
        end do

        if (j == n11) go to 510
        i2 = min0(n11 - j, k1)
        !  transformation to the left hand side with respect to g1.
        do i = 1, i2
          i1 = i + 1
          call fprota(cos, sin, h1(i1), g1(j, i1))
          h1(i) = h1(i1)
        end do

        h1(i1) = 0.
      end do

      !  rotation with the rows n11+1,...n7
510   do j = 1, k1
        ij = n11 + j
        if (ij <= 0) cycle
        piv = h2(j)
        !  calculate the parameters of the givens transformation
        call fpgivs(piv, g2(ij, j), cos, sin)
        !  transformation to the right hand side.
        j1 = ij
        do j2 = 1, idim
          call fprota(cos, sin, xi(j2), c(j1))
          j1 = j1 + n
        end do

        if (j == k1) exit
        j1 = j + 1
        !  transformation to the left hand side.
        do i = j1, k1
          call fprota(cos, sin, h2(i), g2(ij, i))
        end do

      end do

    end do

    !  backward substitution to obtain the b-spline coefficients
    j1 = 1
    do j2 = 1, idim
      call fpbacp(g1, g2, c(j1), n7, k1, c(j1), k2, nest)
      j1 = j1 + n
    end do

    !  calculate from condition (**) the remaining b-spline coefficients.
    do i = 1, k
      j1 = i
      do j = 1, idim
        j2 = j1 + n7
        c(j2) = c(j1)
        j1 = j1 + n
      end do
    end do

    !  computation of f(p).
    fp = 0.
    l = k1
    jj = 0
    do it = 1, m1
      if (u(it) < t(l)) go to 550
      l = l + 1
550   l0 = l - k2
      term = 0.
      do j2 = 1, idim
        fac = 0.
        j1 = l0
        do j = 1, k1
          j1 = j1 + 1
          fac = fac + c(j1) * q(it, j)
        end do

        jj = jj + 1
        term = term + (fac - x(jj))**2
        l0 = l0 + n
      end do

      fp = fp + term * w(it)**2
    end do

    !  test whether the approximation sp(u) is an acceptable solution.
    fpms = fp - s
    if (abs(fpms) < acc) return
    !  test whether the maximal number of iterations is reached.
    if (iter == maxit) then
      ier = 3
      return
    end if

    !  carry out one more step of the iteration process.
    p2 = p
    f2 = fpms
    if (ich3 /= 0) go to 580
    if ((f2 - f3) > acc) go to 575
    !  our initial choice of p is too large.
    p3 = p2
    f3 = f2
    p = p * con4
    if (p <= p1) p = p1 * con9 + p2 * con1
    cycle
575 if (f2 < 0.) ich3 = 1
580 if (ich1 /= 0) go to 590
    if ((f1 - f2) > acc) go to 585
    !  our initial choice of p is too small
    p1 = p2
    f1 = f2
    p = p / con4
    if (p3 < 0._8) cycle
    if (p >= p3) p = p2 * con1 + p3 * con9
    cycle
585 if (f2 > 0._8) ich1 = 1
    !  test whether the iteration process proceeds as theoretically
    !  expected.
590 if (f2 >= f1 .or. f2 <= f3) then
      ier = 2
      return
    end if

    !  find the new value for p.
    p = fprati(p1, f1, p2, f2, p3, f3)
  end do

  !  error codes and messages.
640 ier = -2
  !  the point (z(1),z(2),...,z(idim)) is a solution of our problem.
  !  a constant function is a spline of degree k with all b-spline
  !  coefficients equal to that constant.
  do i = 1, k1
    rn = k1 - i
    t(i) = u(1) - rn * per
    j = i + k1
    rn = i - 1
    t(j) = u(m) + rn * per
  end do

  n = nmin
  j1 = 0
  do j = 1, idim
    fac = z(j)
    j2 = j1
    do i = 1, k1
      j2 = j2 + 1
      c(j2) = fac
    end do

    j1 = j1 + n
  end do

  fp = fp0
  fpint(n) = fp0
  fpint(n - 1) = 0._8
  nrdata(n) = 0

end subroutine fpclos

subroutine fpcurf(iopt, x, y, w, m, xb, xe, k, s, nest, tol, maxit, k1, k2,&
  &n, t, c, fp, fpint, z, a, b, g, q, nrdata, ier)
  !  ..
  !  ..scalar arguments..
  real(8) :: xb, xe, s, tol, fp
  integer :: iopt, m, k, nest, maxit, k1, k2, n, ier
  !  ..array arguments..
  real(8) :: x(m), y(m), w(m), t(nest), c(nest), fpint(nest),&
    &z(nest), a(nest, k1), b(nest, k2), g(nest, k2), q(m, k1)
  integer :: nrdata(nest)
  !  ..local scalars..
  real(8) :: acc, con1, con4, con9, cos, half, fpart, fpms, fpold, fp0, f1, f2, f3
  real(8) :: one, p, pinv, piv, p1, p2, p3, rn, sin, store, term, wi, xi, yi
  integer :: i, ich1, ich3, it, iter, i1, i2, i3, j, k3, l, l0
  integer :: mk1, new, nk1, nmax, nmin, nplus, npl1, nrint, n8
  !  ..local arrays..
  real(8) :: h(7)
  !  ..function references
  ! real(8) :: abs, fprati
  !  ..subroutine references..
  !    fpback,fpbspl,fpgivs,fpdisc,fpknot,fprota
  !  ..
  !  set constants
  one = 1._8
  con1 = 0.1_8
  con9 = 0.9_8
  con4 = 0.04_8
  half = 0.5_8
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  part 1: determination of the number of knots and their position     c
  !  **************************************************************      c
  !  given a set of knots we compute the least-squares spline sinf(x),   c
  !  and the corresponding sum of squared residuals fp=f(p=inf).         c
  !  if iopt=-1 sinf(x) is the requested approximation.                  c
  !  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
  !    if fp <=s we will continue with the current set of knots.         c
  !    if fp > s we will increase the number of knots and compute the    c
  !       corresponding least-squares spline until finally fp<=s.        c
  !    the initial choice of knots depends on the value of s and iopt.   c
  !    if s=0 we have spline interpolation; in that case the number of   c
  !    knots equals nmax = m+k+1.                                        c
  !    if s > 0 and                                                      c
  !      iopt=0 we first compute the least-squares polynomial of         c
  !      degree k; n = nmin = 2*k+2                                      c
  !      iopt=1 we start with the set of knots found at the last         c
  !      call of the routine, except for the case that s > fp0; then     c
  !      we compute directly the least-squares polynomial of degree k.   c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  determine nmin, the number of knots for polynomial approximation.
  ! JF: Initialize coefficients
  ! JF: It seems that this is needed!!
  c = 0._8

  nmin = 2 * k1
  if (iopt < 0) go to 60
  !  calculation of acc, the absolute tolerance for the root of f(p)=s.
  acc = tol * s
  !  determine nmax, the number of knots for spline interpolation.
  nmax = m + k1
  if (s > 0._8) go to 45
  !  if s=0, s(x) is an interpolating spline.
  !  test whether the required storage space exceeds the available one.
  n = nmax
  if (nmax > nest) go to 420
  !  find the position of the interior knots in case of interpolation.
10 mk1 = m - k1
  if (mk1 == 0) go to 60
  k3 = k / 2
  i = k2
  j = k3 + 2
  if (k3 * 2 == k) go to 30
  do l = 1, mk1
    t(i) = x(j)
    i = i + 1
    j = j + 1
  end do

  go to 60
30 do l = 1, mk1
    t(i) = (x(j) + x(j - 1)) * half
    i = i + 1
    j = j + 1
  end do
  go to 60
  !  if s>0 our initial choice of knots depends on the value of iopt.
  !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
  !  polynomial of degree k which is a spline without interior knots.
  !  if iopt=1 and fp0>s we start computing the least squares spline
  !  according to the set of knots found at the last call of the routine.
45 if (iopt == 0) go to 50
  if (n == nmin) go to 50
  fp0 = fpint(n)
  fpold = fpint(n - 1)
  nplus = nrdata(n)
  if (fp0 > s) go to 60
50 n = nmin
  fpold = 0._8
  nplus = 0
  nrdata(1) = m - 2

  !JF: Not needed?
  ! a = 0._8

  !  main loop for the different sets of knots. m is a save upper bound
  !  for the number of trials.
60 do iter = 1, m
    if (n == nmin) ier = -2
    !  find nrint, tne number of knot intervals.
    nrint = n - nmin + 1
    !  find the position of the additional knots which are needed for
    !  the b-spline representation of s(x).
    nk1 = n - k1
    i = n
    do j = 1, k1
      t(j) = xb
      t(i) = xe
      i = i - 1
    end do

    !  compute the b-spline coefficients of the least-squares spline
    !  sinf(x). the observation matrix a is built up row by row and
    !  reduced to upper triangular form by givens transformations.
    !  at the same time fp=f(p=inf) is computed.
    fp = 0._8
    !  initialize the observation matrix a.

    z(:nk1) = 0._8
    a(:nk1, :k1) = 0._8
    l = k1
    do it = 1, m
      !  fetch the current data point x(it),y(it).
      xi = x(it)
      wi = w(it)
      yi = y(it) * wi
      !  search for knot interval t(l) <= xi < t(l+1).
      do while (t(l + 1) <= xi .and. l /= nk1)
        l = l + 1
      end do
      !  evaluate the (k+1) non-zero b-splines at xi and store them in q.
      call fpbspl(t, n, k, xi, l, h)
      q(it, :k1) = h(:k1)
      h(:k1) = h(:k1) * wi
      !  rotate the new row of the observation matrix into triangle.
      j = l - k1
      do i = 1, k1
        j = j + 1
        piv = h(i)
        if (piv == 0._8) cycle
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, a(j, 1), cos, sin)
        !  transformations to right hand side.
        call fprota(cos, sin, yi, z(j))
        if (i == k1) exit
        i2 = 1
        i3 = i + 1
        do i1 = i3, k1
          i2 = i2 + 1
          !  transformations to left hand side.
          call fprota(cos, sin, h(i1), a(j, i2))
        end do
      end do

      !  add contribution of this row to the sum of squares of residual
      !  right hand sides.
      fp = fp + yi * yi
    end do

    if (ier == (-2)) fp0 = fp
    fpint(n) = fp0
    fpint(n - 1) = fpold
    nrdata(n) = nplus
    !  backward substitution to obtain the b-spline coefficients.
    call fpback(a, z, nk1, k1, c, nest)
    !  test whether the approximation sinf(x) is an acceptable solution.
    if (iopt < 0) go to 440
    fpms = fp - s
    if (abs(fpms) < acc) go to 440
    !  if f(p=inf) < s accept the choice of knots.
    if (fpms < 0._8) go to 250
    !  if n = nmax, sinf(x) is an interpolating spline.
    if (n == nmax) go to 430
    !  increase the number of knots.
    !  if n=nest we cannot increase the number of knots because of
    !  the storage capacity limitation.
    if (n == nest) go to 420
    !  determine the number of knots nplus we are going to add.
    if (ier == 0) go to 140
    nplus = 1
    ier = 0
    go to 150
140 npl1 = nplus * 2
    rn = nplus
    if (fpold - fp > acc) npl1 = int(rn * fpms / (fpold - fp))
    nplus = min(nplus * 2, max(npl1, nplus / 2, 1))
150 fpold = fp
    !  compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
    !  t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
    fpart = 0._8
    i = 1
    l = k2
    new = 0
    do it = 1, m
      if (x(it) >= t(l) .and. l <= nk1) then
        new = 1
        l = l + 1
      end if

      term = 0._8
      l0 = l - k2
      do j = 1, k1
        l0 = l0 + 1
        term = term + c(l0) * q(it, j)
      end do

      term = (w(it) * (term - y(it)))**2
      fpart = fpart + term
      if (new == 0) cycle
      store = term * half
      fpint(i) = fpart - store
      i = i + 1
      fpart = store
      new = 0
    end do

    fpint(nrint) = fpart
    do l = 1, nplus
      !  add a new knot.
      call fpknot(x, m, t, n, fpint, nrdata, nrint, nest, 1)
      !  if n=nmax we locate the knots as for interpolation.
      if (n == nmax) go to 10
      !  test whether we cannot further increase the number of knots.
      if (n == nest) exit
    end do

    !  restart the computations with the new set of knots.
  end do

  !  test whether the least-squares kth degree polynomial is a solution
  !  of our approximation problem.
250 if (ier == (-2)) go to 440
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  part 2: determination of the smoothing spline sp(x).                c
  !  ***************************************************                 c
  !  we have determined the number of knots and their position.          c
  !  we now compute the b-spline coefficients of the smoothing spline    c
  !  sp(x). the observation matrix a is extended by the rows of matrix   c
  !  b expressing that the kth derivative discontinuities of sp(x) at    c
  !  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
  !  ponding weights of these additional rows are set to 1/p.            c
  !  iteratively we then have to determine the value of p such that      c
  !  f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that    c
  !  the least-squares kth degree polynomial corresponds to p=0, and     c
  !  that the least-squares spline corresponds to p=infinity. the        c
  !  iteration process which is proposed here, makes use of rational     c
  !  interpolation. since f(p) is a convex and strictly decreasing       c
  !  function of p, it can be approximated by a rational function        c
  !  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
  !  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
  !  to calculate the new value of p such that r(p)=s. convergence is    c
  !  guaranteed by taking f1>0 and f3<0.                                 c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  evaluate the discontinuity jump of the kth derivative of the
  !  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
  call fpdisc(t, n, k2, b, nest)
  !  initial value for p.
  p1 = 0._8
  f1 = fp0 - s
  p3 = -one
  f3 = fpms
  p = sum(a(:nk1, 1))

  rn = nk1
  p = rn / p
  ich1 = 0
  ich3 = 0
  n8 = n - nmin
  !  iteration process to find the root of f(p) = s.
  do iter = 1, maxit
    !  the rows of matrix b with weight 1/p are rotated into the
    !  triangularised observation matrix a which is stored in g.
    pinv = one / p
    c(:nk1) = z(:nk1)
    g(:nk1, k2) = 0._8
    g(:nk1, :k1) = a(:nk1, :k1)
    do it = 1, n8
      !  the row of matrix b is rotated into triangle by givens transformation
      h(:k2) = b(it, :k2) * pinv
      yi = 0._8
      do j = it, nk1
        piv = h(1)
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, g(j, 1), cos, sin)
        !  transformations to right hand side.
        call fprota(cos, sin, yi, c(j))
        if (j == nk1) exit
        i2 = k1; IF (j > n8) i2 = nk1 - j
        do i = 1, i2
          !  transformations to left hand side.
          i1 = i + 1
          call fprota(cos, sin, h(i1), g(j, i1))
          h(i) = h(i1)
        end do
        h(i2 + 1) = 0._8
      end do
    end do
    !  backward substitution to obtain the b-spline coefficients.
    call fpback(g, c, nk1, k2, c, nest)
    !  computation of f(p).
    fp = 0._8
    l = k2
    do it = 1, m
      if (x(it) >= t(l) .and. l <= nk1) l = l + 1
      l0 = l - k2
      term = 0._8
      do j = 1, k1
        l0 = l0 + 1
        term = term + c(l0) * q(it, j)
      end do

      fp = fp + (w(it) * (term - y(it)))**2
    end do

    !  test whether the approximation sp(x) is an acceptable solution.
    fpms = fp - s
    if (abs(fpms) < acc) go to 440
    !  test whether the maximal number of iterations is reached.
    if (iter == maxit) go to 400
    !  carry out one more step of the iteration process.
    p2 = p
    f2 = fpms
    if (ich3 /= 0) go to 340
    if ((f2 - f3) > acc) go to 335
    !  our initial choice of p is too large.
    p3 = p2
    f3 = f2
    p = p * con4
    if (p <= p1) p = p1 * con9 + p2 * con1
    cycle
335 if (f2 < 0._8) ich3 = 1
340 if (ich1 /= 0) go to 350
    if ((f1 - f2) > acc) go to 345
    !  our initial choice of p is too small
    p1 = p2
    f1 = f2
    p = p / con4
    if (p3 < 0.) cycle
    if (p >= p3) p = p2 * con1 + p3 * con9
    cycle
345 if (f2 > 0._8) ich1 = 1
    !  test whether the iteration process proceeds as theoretically
    !  expected.
350 if (f2 >= f1 .or. f2 <= f3) go to 410
    !  find the new value for p.
    p = fprati(p1, f1, p2, f2, p3, f3)
  end do

  !  error codes and messages.
400 ier = 3
  go to 440
410 ier = 2
  go to 440
420 ier = 1
  go to 440
430 ier = -1
440 return
end subroutine fpcurf

subroutine fpdisc(t, n, k2, b, nest)
  !  subroutine fpdisc calculates the discontinuity jumps of the kth
  !  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
  !  ..scalar arguments..
  integer :: n, k2, nest
  !  ..array arguments..
  real(8) :: t(n), b(nest, k2)
  !  ..local scalars..
  real(8) :: an, fac, prod
  integer :: i, ik, j, jk, k, k1, l, lj, lk, lmk, lp, nk1, nrint
  !  ..local array..
  real(8) :: h(12)
  !  ..
  k1 = k2 - 1
  k = k1 - 1
  nk1 = n - k1
  nrint = nk1 - k
  an = nrint
  fac = an / (t(nk1 + 1) - t(k1))
  do l = k2, nk1
    lmk = l - k1
    do j = 1, k1
      ik = j + k1
      lj = l + j
      lk = lj - k2
      h(j) = t(l) - t(lk)
      h(ik) = t(l) - t(lj)
    end do

    lp = lmk
    do j = 1, k2
      jk = j
      prod = h(j)
      do i = 1, k
        jk = jk + 1
        prod = prod * h(jk) * fac
      end do

      lk = lp + k1
      b(lmk, j) = (t(lk) - t(lp)) / prod
      lp = lp + 1
    end do

  end do
end subroutine fpdisc

subroutine fpgivs(piv, ww, cos, sin)
  !  subroutine fpgivs calculates the parameters of a givens
  !  transformation .
  !  ..
  !  ..scalar arguments..
  real(8) :: piv, ww, cos, sin
  !  ..local scalars..
  real(8) :: dd, one, store
  !  ..function references..
  ! real(8) :: abs, sqrt
  !  ..
  one = 1._8
  store = abs(piv)
  if (store >= ww) dd = store * sqrt(one + (ww / piv)**2)
  if (store < ww) dd = ww * sqrt(one + (piv / ww)**2)
  cos = ww / dd
  sin = piv / dd
  ww = dd
end subroutine fpgivs

subroutine fpknot(x, m, t, n, fpint, nrdata, nrint, nest, istart)
  implicit none
  !  subroutine fpknot locates an additional knot for a spline of degree
  !  k and adjusts the corresponding parameters,i.e.
  !    t     : the position of the knots.
  !    n     : the number of knots.
  !    nrint : the number of knotintervals.
  !    fpint : the sum of squares of residual right hand sides
  !            for each knot interval.
  !    nrdata: the number of data points inside each knot interval.
  !  istart indicates that the smallest data point at which the new knot
  !  may be added is x(istart+1)
  !  ..
  !  ..scalar arguments..
  integer :: m, n, nrint, nest, istart
  !  ..array arguments..
  real(8) :: x(m), t(nest), fpint(nest)
  integer :: nrdata(nest)
  !  ..local scalars..
  real(8) :: an, am, fpmax
  integer :: ihalf, j, jbegin, jj, jk, jpoint, k, maxbeg, maxpt, next, nrx, number
  !  ..
  k = (n - nrint - 1) / 2
  !  search for knot interval t(number+k) <= x <= t(number+k+1) where
  !  fpint(number) is maximal on the condition that nrdata(number)
  !  not equals zero.
  fpmax = 0.
  jbegin = istart
  do j = 1, nrint
    jpoint = nrdata(j)
    if (fpmax < fpint(j) .and. jpoint /= 0) then
      fpmax = fpint(j)
      number = j
      maxpt = jpoint
      maxbeg = jbegin
    end if
    jbegin = jbegin + jpoint + 1
  end do
  !  let coincide the new knot t(number+k+1) with a data point x(nrx)
  !  inside the old knot interval t(number+k) <= x <= t(number+k+1).
  ihalf = maxpt / 2 + 1
  nrx = maxbeg + ihalf
  next = number + 1
  if (next <= nrint) then
    !  adjust the different parameters.
    do j = next, nrint
      jj = next + nrint - j
      fpint(jj + 1) = fpint(jj)
      nrdata(jj + 1) = nrdata(jj)
      jk = jj + k
      t(jk + 1) = t(jk)
    end do
  end if

  nrdata(number) = ihalf - 1
  nrdata(next) = maxpt - ihalf
  am = maxpt
  an = nrdata(number)
  fpint(number) = fpmax * an / am
  an = nrdata(next)
  fpint(next) = fpmax * an / am
  jk = next + k
  t(jk) = x(nrx)
  n = n + 1
  nrint = nrint + 1

end subroutine fpknot

subroutine fppara(iopt, idim, m, u, mx, x, w, ub, ue, k, s, nest, tol, maxit,&
  &k1, k2, n, t, nc, c, fp, fpint, z, a, b, g, q, nrdata, ier)
  !  ..
  !  ..scalar arguments..
  real(8) :: ub, ue, s, tol, fp
  integer :: iopt, idim, m, mx, k, nest, maxit, k1, k2, n, nc, ier
  !  ..array arguments..
  real(8) :: u(m), x(mx), w(m), t(nest), c(nc), fpint(nest),&
    &z(nc), a(nest, k1), b(nest, k2), g(nest, k2), q(m, k1)
  integer :: nrdata(nest)
  !  ..local scalars..
  real(8) :: acc, con1, con4, con9, cos, fac, fpart, fpms, fpold, fp0, f1, f2, f3,&
    &half, one, p, pinv, piv, p1, p2, p3, rn, sin, store, term, ui, wi
  integer :: i, ich1, ich3, it, iter, i1, i2, i3, j, jj, j1, j2, k3, l, l0,&
    &mk1, new, nk1, nmax, nmin, nplus, npl1, nrint, n8
  !  ..local arrays..
  real(8) :: h(7), xi(10)
  !  ..function references
  ! real(8) :: abs, fprati
  !  ..subroutine references..
  !    fpback,fpbspl,fpgivs,fpdisc,fpknot,fprota
  !  ..
  !  set constants
  one = 1._8
  con1 = 0.1_8
  con9 = 0.9_8
  con4 = 0.04_8
  half = 0.5_8
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  part 1: determination of the number of knots and their position     c
  !  **************************************************************      c
  !  given a set of knots we compute the least-squares curve sinf(u),    c
  !  and the corresponding sum of squared residuals fp=f(p=inf).         c
  !  if iopt=-1 sinf(u) is the requested curve.                          c
  !  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
  !    if fp <=s we will continue with the current set of knots.         c
  !    if fp > s we will increase the number of knots and compute the    c
  !       corresponding least-squares curve until finally fp<=s.         c
  !    the initial choice of knots depends on the value of s and iopt.   c
  !    if s=0 we have spline interpolation; in that case the number of   c
  !    knots equals nmax = m+k+1.                                        c
  !    if s > 0 and                                                      c
  !      iopt=0 we first compute the least-squares polynomial curve of   c
  !      degree k; n = nmin = 2*k+2                                      c
  !      iopt=1 we start with the set of knots found at the last         c
  !      call of the routine, except for the case that s > fp0; then     c
  !      we compute directly the polynomial curve of degree k.           c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  determine nmin, the number of knots for polynomial approximation.
  nmin = 2 * k1
  if (iopt < 0) go to 60
  !  calculation of acc, the absolute tolerance for the root of f(p)=s.
  acc = tol * s
  !  determine nmax, the number of knots for spline interpolation.
  nmax = m + k1
  if (s > 0._8) go to 45        !  S(u) is an interpolating curve.

  n = nmax
  if (nmax > nest) then !  The required storage space exceeds the available one?
    ier = 1
    return
  end if

  !  find the position of the interior knots in case of interpolation.
10 mk1 = m - k1
  if (mk1 /= 0) then
    k3 = k / 2
    i = k2
    j = k3 + 2
    if (k3 * 2 == k) then   ! Degree k is even. Take the center points
      t(k2:k2 + mk1 - 1) = (u(k3 + 2:k3 + 1 + mk1) + u(k3 + 1:k3 + mk1)) * half
    else                    ! Odd values of k
      t(k2:k2 + mk1 - 1) = u(k3 + 2:k3 + 1 + mk1)
    end if
  end if
  go to 60

  !  if s>0 our initial choice of knots depends on the value of iopt.
  !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
  !  polynomial curve which is a spline curve without interior knots.
  !  if iopt=1 and fp0>s we start computing the least squares spline curve
  !  according to the set of knots found at the last call of the routine.
45 if (iopt == 0) go to 50
  if (n == nmin) go to 50
  fp0 = fpint(n)
  fpold = fpint(n - 1)
  nplus = nrdata(n)
  if (fp0 > s) go to 60
50 n = nmin
  fpold = 0._8
  nplus = 0
  nrdata(1) = m - 2
  !  main loop for the different sets of knots. m is a safe upper bound
  !  for the number of trials.
60 do iter = 1, m
    if (n == nmin) ier = -2
    !  find nrint, tne number of knot intervals.
    nrint = n - nmin + 1
    !  find the position of the additional knots which are needed for
    !  the b-spline representation of s(u).
    nk1 = n - k1
    t(1:k1) = ub
    t(n - k:n) = ue
    !  compute the b-spline coefficients of the least-squares spline curve
    !  sinf(u). the observation matrix a is built up row by row and
    !  reduced to upper triangular form by givens transformations.
    !  at the same time fp=f(p=inf) is computed.
    fp = 0._8
    !  initialize the b-spline coefficients and the observation matrix a.
    z(:nc) = 0._8
    a(:nk1, :k1) = 0._8
    l = k1
    do it = 1, m
      !  fetch the current data point u(it),x(it).
      ui = u(it)
      wi = w(it)
      jj = (it - 1) * idim             ! Temporary value
      xi(1:idim) = x(jj + 1:jj + idim) * wi
      !  search for knot interval t(l) <= ui < t(l+1).
      do while (ui >= t(l + 1) .and. l < nk1)
        l = l + 1
      end do
      !  evaluate the (k+1) non-zero b-splines at ui and store them in q.
      call fpbspl(t, n, k, ui, l, h)
      q(it, :k1) = h(:k1)
      h(:k1) = h(:k1) * wi
      !  rotate the new row of the observation matrix into triangle.
      j = l - k1
      do i = 1, k1
        j = j + 1
        piv = h(i)
        if (piv == 0._8) cycle
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, a(j, 1), cos, sin)
        !  transformations to right hand side.
        j1 = j
        do j2 = 1, idim
          call fprota(cos, sin, xi(j2), z(j1))
          j1 = j1 + n
        end do

        if (i == k1) exit
        i2 = 1
        i3 = i + 1
        do i1 = i3, k1
          i2 = i2 + 1
          !  transformations to left hand side.
          call fprota(cos, sin, h(i1), a(j, i2))
        end do
      end do
      !  add contribution of this row to the sum of squares of residual
      !  right hand sides.
      do j2 = 1, idim
        fp = fp + xi(j2)**2
      end do
    end do

    if (ier == (-2)) fp0 = fp
    fpint(n) = fp0
    fpint(n - 1) = fpold
    nrdata(n) = nplus
    !  backward substitution to obtain the b-spline coefficients.
    j1 = 1
    do j2 = 1, idim
      call fpback(a, z(j1), nk1, k1, c(j1), nest)
      j1 = j1 + n
    end do

    !  test whether the approximation sinf(u) is an acceptable solution.
    if (iopt < 0) return
    fpms = fp - s
    if (abs(fpms) < acc) return
    !  if f(p=inf) < s accept the choice of knots.
    if (fpms < 0._8) go to 250
    !  if n = nmax, sinf(u) is an interpolating spline curve.
    if (n == nmax) then
      ier = -1
      return
    end if
    !  increase the number of knots.
    !  if n=nest we cannot increase the number of knots because of
    !  the storage capacity limitation.
    if (n == nest) then
      ier = 1
      return
    end if
    !  determine the number of knots nplus we are going to add.
    if (ier == 0) then
      npl1 = nplus * 2
      rn = nplus
      if (fpold - fp > acc) npl1 = int(rn * fpms / (fpold - fp))
      nplus = min(nplus * 2, max(npl1, nplus / 2, 1))
    else
      nplus = 1
      ier = 0
    end if
    fpold = fp
    !  compute the sum of squared residuals for each knot interval
    !  t(j+k) <= u(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
    fpart = 0._8
    i = 1
    l = k2
    new = 0
    jj = 0
    do it = 1, m
      if (u(it) >= t(l) .and. l <= nk1) then
        new = 1
        l = l + 1
      end if
      term = 0._8
      l0 = l - k2
      do j2 = 1, idim
        fac = sum(c(l0 + 1:l0 + k1) * q(it, 1:k1))
        jj = jj + 1
        term = term + (w(it) * (fac - x(jj)))**2
        l0 = l0 + n
      end do

      fpart = fpart + term
      if (new /= 0) then
        store = term * half
        fpint(i) = fpart - store
        i = i + 1
        fpart = store
        new = 0
      end if
    end do

    fpint(nrint) = fpart
    do l = 1, nplus
      !  add a new knot.
      call fpknot(u, m, t, n, fpint, nrdata, nrint, nest, 1)
      !  if n=nmax we locate the knots as for interpolation
      if (n == nmax) go to 10
      !  test whether we cannot further increase the number of knots.
      if (n == nest) exit
    end do
    !  restart the computations with the new set of knots.
  end do
  !  test whether the least-squares kth degree polynomial curve is a
  !  solution of our approximation problem.
250 if (ier == (-2)) return
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  part 2: determination of the smoothing spline curve sp(u).          c
  !  **********************************************************          c
  !  we have determined the number of knots and their position.          c
  !  we now compute the b-spline coefficients of the smoothing curve     c
  !  sp(u). the observation matrix a is extended by the rows of matrix   c
  !  b expressing that the kth derivative discontinuities of sp(u) at    c
  !  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
  !  ponding weights of these additional rows are set to 1/p.            c
  !  iteratively we then have to determine the value of p such that f(p),c
  !  the sum of squared residuals be = s. we already know that the least c
  !  squares kth degree polynomial curve corresponds to p=0, and that    c
  !  the least-squares spline curve corresponds to p=infinity. the       c
  !  iteration process which is proposed here, makes use of rational     c
  !  interpolation. since f(p) is a convex and strictly decreasing       c
  !  function of p, it can be approximated by a rational function        c
  !  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
  !  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
  !  to calculate the new value of p such that r(p)=s. convergence is    c
  !  guaranteed by taking f1>0 and f3<0.                                 c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  evaluate the discontinuity jump of the kth derivative of the
  !  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
  call fpdisc(t, n, k2, b, nest)
  !  initial value for p.
  p1 = 0.
  f1 = fp0 - s
  p3 = -one
  f3 = fpms
  p = sum(a(:nk1, 1))
  rn = nk1
  p = rn / p
  ich1 = 0
  ich3 = 0
  n8 = n - nmin
  !  iteration process to find the root of f(p) = s.
  do iter = 1, maxit
    !  the rows of matrix b with weight 1/p are rotated into the
    !  triangularised observation matrix a which is stored in g.
    pinv = one / p
    c(:nc) = z(:nc)
    g(:nk1, k2) = 0._8
    g(:nk1, :k1) = a(:nk1, :k1)

    do it = 1, n8
      !  the row of matrix b is rotated into triangle by givens transformation
      h(:k2) = b(it, :k2) * pinv
      xi(:idim) = 0._8

      inner1: do j = it, nk1
        piv = h(1)
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, g(j, 1), cos, sin)
        !  transformations to right hand side.
        j1 = j
        do j2 = 1, idim
          call fprota(cos, sin, xi(j2), c(j1))
          j1 = j1 + n
        end do

        if (j == nk1) exit inner1
        i2 = k1
        if (j > n8) i2 = nk1 - j
        do i = 1, i2
          !  transformations to left hand side.
          i1 = i + 1
          call fprota(cos, sin, h(i1), g(j, i1))
          h(i) = h(i1)
        end do

        h(i2 + 1) = 0.
      end do inner1
    end do

    !  backward substitution to obtain the b-spline coefficients.
    j1 = 1
    do j2 = 1, idim
      call fpback(g, c(j1), nk1, k2, c(j1), nest)
      j1 = j1 + n
    end do

    !  computation of f(p).
    fp = 0.
    l = k2
    jj = 0
    do it = 1, m
      if ((u(it) >= t(l)) .and. (l <= nk1)) l = l + 1
      l0 = l - k2
      term = 0.
      do j2 = 1, idim
        fac = sum(c(l0 + 1:l0 + k1) * q(it, 1:k1))

        jj = jj + 1
        term = term + (fac - x(jj))**2
        l0 = l0 + n
      end do

      fp = fp + term * w(it)**2
    end do
    !  test whether the approximation sp(u) is an acceptable solution.
    fpms = fp - s
    if (abs(fpms) < acc) return
    !  test whether the maximal number of iterations is reached.
    if (iter == maxit) then
      ier = 3
      return
    end if

    !  carry out one more step of the iteration process.
    p2 = p
    f2 = fpms
    if (ich3 == 0) then
      if ((f2 - f3) <= acc) then
        !  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p * con4
        if (p <= p1) p = p1 * con9 + p2 * con1
        cycle
      end if

      IF (f2 < 0._8) ich3 = 1
    end if

    if (ich1 == 0) then
      if ((f1 - f2) <= acc) then
        !  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p / con4
        if (p3 >= 0._8) then
          IF (p >= p3) p = p2 * con1 + p3 * con9
        end if
        cycle
      end if

      IF (f2 > 0._8) ich1 = 1
    end if

    !  Test whether the iteration process proceeds as expected.
    if (f2 >= f1 .or. f2 <= f3) then
      ier = 2
      return
    end if

    !  find the new value for p.
    p = fprati(p1, f1, p2, f2, p3, f3)
  end do
end subroutine fppara

subroutine fpperi(iopt, x, y, w, m, k, s, nest, tol, maxit, k1, k2, n, t, c,&
  &fp, fpint, z, a1, a2, b, g1, g2, q, nrdata, ier)
  !  ..
  !  ..scalar arguments..
  real(8) :: s, tol, fp
  integer :: iopt, m, k, nest, maxit, k1, k2, n, ier
  !  ..array arguments..
  real(8) :: x(m), y(m), w(m), t(nest), c(nest), fpint(nest), z(nest),&
    &a1(nest, k1), a2(nest, k), b(nest, k2), g1(nest, k2), g2(nest, k1),&
    &q(m, k1)
  integer :: nrdata(nest)
  !  ..local scalars..
  real(8) :: acc, cos, c1, d1, fpart, fpms, fpold, fp0, f1, f2, f3, p, per, pinv, piv,&
    &p1, p2, p3, sin, store, term, wi, xi, yi, rn, one, con1, con4, con9, half
  integer :: i, ich1, ich3, ij, ik, it, iter, i1, i2, i3, j, jk, jper, j1, j2, kk,&
    &kk1, k3, l, l0, l1, l5, mm, m1, new, nk1, nk2, nmax, nmin, nplus, npl1,&
    &nrint, n10, n11, n7, n8
  !  ..local arrays..
  real(8) :: h(6), h1(7), h2(6)
  !  ..function references..
  ! real(8) :: abs, fprati
  integer :: max0, min0
  !  ..subroutine references..
  !    fpbacp,fpbspl,fpgivs,fpdisc,fpknot,fprota
  !  ..
  !  set constants
  one = 0.1e+01
  con1 = 0.1e0
  con9 = 0.9e0
  con4 = 0.4e-01
  half = 0.5e0
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  part 1: determination of the number of knots and their position     c
  !  **************************************************************      c
  !  given a set of knots we compute the least-squares periodic spline   c
  !  sinf(x). if the sum f(p=inf) <= s we accept the choice of knots.    c
  !  the initial choice of knots depends on the value of s and iopt.     c
  !    if s=0 we have spline interpolation; in that case the number of   c
  !    knots equals nmax = m+2*k.                                        c
  !    if s > 0 and                                                      c
  !      iopt=0 we first compute the least-squares polynomial of         c
  !      degree k; n = nmin = 2*k+2. since s(x) must be periodic we      c
  !      find that s(x) is a constant function.                          c
  !      iopt=1 we start with the set of knots found at the last         c
  !      call of the routine, except for the case that s > fp0; then     c
  !      we compute directly the least-squares periodic polynomial.      c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  m1 = m - 1
  kk = k
  kk1 = k1
  k3 = 3 * k + 1
  nmin = 2 * k1
  !  determine the length of the period of s(x).
  per = x(m) - x(1)
  if (iopt < 0) go to 50
  !  calculation of acc, the absolute tolerance for the root of f(p)=s.
  acc = tol * s
  !  determine nmax, the number of knots for periodic spline interpolation
  nmax = m + 2 * k
  if (s > 0. .or. nmax == nmin) go to 30
  !  if s=0, s(x) is an interpolating spline.
  n = nmax
  !  test whether the required storage space exceeds the available one.
  if (n > nest) go to 620
  !  find the position of the interior knots in case of interpolation.
5 if ((k / 2) * 2 == k) go to 20
  do i = 2, m1
    j = i + k
    t(j) = x(i)
  end do

  if (s > 0.) go to 50
  kk = k - 1
  kk1 = k
  if (kk > 0) go to 50
  t(1) = t(m) - per
  t(2) = x(1)
  t(m + 1) = x(m)
  t(m + 2) = t(3) + per
  c(:m1) = y(:m1)
  c(m) = c(1)
  fp = 0.
  fpint(n) = fp0
  fpint(n - 1) = 0.
  nrdata(n) = 0
  go to 630
20 do i = 2, m1
    j = i + k
    t(j) = (x(i) + x(i - 1)) * half
  end do

  go to 50
  !  if s > 0 our initial choice depends on the value of iopt.
  !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
  !  periodic polynomial. (i.e. a constant function).
  !  if iopt=1 and fp0>s we start computing the least-squares periodic
  !  spline according the set of knots found at the last call of the
  !  routine.
30 if (iopt == 0) go to 35
  if (n == nmin) go to 35
  fp0 = fpint(n)
  fpold = fpint(n - 1)
  nplus = nrdata(n)
  if (fp0 > s) go to 50
  !  the case that s(x) is a constant function is treated separetely.
  !  find the least-squares constant c1 and compute fp0 at the same time.
35 fp0 = 0.
  d1 = 0.
  c1 = 0.
  do it = 1, m1
    wi = w(it)
    yi = y(it) * wi
    call fpgivs(wi, d1, cos, sin)
    call fprota(cos, sin, yi, c1)
    fp0 = fp0 + yi**2
  end do

  c1 = c1 / d1
  !  test whether that constant function is a solution of our problem.
  fpms = fp0 - s
  if (fpms < acc .or. nmax == nmin) go to 640
  fpold = fp0
  !  test whether the required storage space exceeds the available one.
  if (nmin >= nest) go to 620
  !  start computing the least-squares periodic spline with one
  !  interior knot.
  nplus = 1
  n = nmin + 1
  mm = (m + 1) / 2
  t(k2) = x(mm)
  nrdata(1) = mm - 2
  nrdata(2) = m1 - mm
  !  main loop for the different sets of knots. m is a save upper
  !  bound for the number of trials.
50 do iter = 1, m
    !  find nrint, the number of knot intervals.
    nrint = n - nmin + 1
    !  find the position of the additional knots which are needed for
    !  the b-spline representation of s(x). if we take
    !      t(k+1) = x(1), t(n-k) = x(m)
    !      t(k+1-j) = t(n-k-j) - per, j=1,2,...k
    !      t(n-k+j) = t(k+1+j) + per, j=1,2,...k
    !  then s(x) is a periodic spline with period per if the b-spline
    !  coefficients satisfy the following conditions
    !      c(n7+j) = c(j), j=1,...k   (**)   with n7=n-2*k-1.
    t(k1) = x(1)
    nk1 = n - k1
    nk2 = nk1 + 1
    t(nk2) = x(m)
    do j = 1, k
      i1 = nk2 + j
      i2 = nk2 - j
      j1 = k1 + j
      j2 = k1 - j
      t(i1) = t(j1) + per
      t(j2) = t(i2) - per
    end do

    !  compute the b-spline coefficients c(j),j=1,...n7 of the least-squares
    !  periodic spline sinf(x). the observation matrix a is built up row
    !  by row while taking into account condition (**) and is reduced to
    !  triangular form by givens transformations .
    !  at the same time fp=f(p=inf) is computed.
    !  the n7 x n7 triangularised upper matrix a has the form
    !            ! a1 '    !
    !        a = !    ' a2 !
    !            ! 0  '    !
    !  with a2 a n7 x k matrix and a1 a n10 x n10 upper triangular
    !  matrix of bandwidth k+1 ( n10 = n7-k).
    !  initialization.
    z(:nk1) = 0._8
    a1(:nk1, :kk1) = 0._8
    n7 = nk1 - k
    n10 = n7 - kk
    jper = 0
    fp = 0.
    l = k1
    do it = 1, m1
      !  fetch the current data point x(it),y(it)
      xi = x(it)
      wi = w(it)
      yi = y(it) * wi
      !  search for knot interval t(l) <= xi < t(l+1).
      do while (xi >= t(l + 1))
        l = l + 1
      end do

      !  evaluate the (k+1) non-zero b-splines at xi and store them in q.
      call fpbspl(t, n, k, xi, l, h)
      do i = 1, k1
        q(it, i) = h(i)
        h(i) = h(i) * wi
      end do

      l5 = l - k1
      !  test whether the b-splines nj,k+1(x),j=1+n7,...nk1 are all zero at xi
      if (l5 < n10) go to 285
      if (jper /= 0) go to 160
      !  initialize the matrix a2.

      a2(:n7, :kk) = 0._8

      jk = n10 + 1
      do i = 1, kk
        ik = jk
        do j = 1, kk1
          if (ik <= 0) exit
          a2(ik, i) = a1(ik, j)
          ik = ik - 1
        end do

        jk = jk + 1
      end do

      jper = 1
      !  if one of the b-splines nj,k+1(x),j=n7+1,...nk1 is not zero at xi
      !  we take account of condition (**) for setting up the new row
      !  of the observation matrix a. this row is stored in the arrays h1
      !  (the part with respect to a1) and h2 (the part with
      !  respect to a2).
160   h1(:kk1) = 0._8
      h2(:kk) = 0._8

      j = l5 - n10
      do i = 1, kk1
        j = j + 1
        l0 = j
180     l1 = l0 - kk
        if (l1 <= 0) go to 200
        if (l1 <= n10) go to 190
        l0 = l1 - n10
        go to 180
190     h1(l1) = h(i)
        cycle
200     h2(l0) = h2(l0) + h(i)
      end do

      !  rotate the new row of the observation matrix into triangle
      !  by givens transformations.
      if (n10 <= 0) go to 250
      !  rotation with the rows 1,2,...n10 of matrix a.
      do j = 1, n10
        piv = h1(1)
        if (piv /= 0.) go to 214
        do i = 1, kk
          h1(i) = h1(i + 1)
        end do

        h1(kk1) = 0._8
        cycle
        !  calculate the parameters of the givens transformation.
214     call fpgivs(piv, a1(j, 1), cos, sin)
        !  transformation to the right hand side.
        call fprota(cos, sin, yi, z(j))
        !  transformations to the left hand side with respect to a2.
        do i = 1, kk
          call fprota(cos, sin, h2(i), a2(j, i))
        end do

        if (j == n10) go to 250
        i2 = min0(n10 - j, kk)
        !  transformations to the left hand side with respect to a1.
        do i = 1, i2
          i1 = i + 1
          call fprota(cos, sin, h1(i1), a1(j, i1))
          h1(i) = h1(i1)
        end do

        h1(i1) = 0.
      end do

      !  rotation with the rows n10+1,...n7 of matrix a.
250   do j = 1, kk
        ij = n10 + j
        if (ij <= 0) cycle
        piv = h2(j)
        if (piv == 0.) cycle
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, a2(ij, j), cos, sin)
        !  transformations to right hand side.
        call fprota(cos, sin, yi, z(ij))
        if (j == kk) go to 280
        j1 = j + 1
        !  transformations to left hand side.
        do i = j1, kk
          call fprota(cos, sin, h2(i), a2(ij, i))
        end do

      end do

      !  add contribution of this row to the sum of squares of residual
      !  right hand sides.
280   fp = fp + yi**2
      cycle
      !  rotation of the new row of the observation matrix into
      !  triangle in case the b-splines nj,k+1(x),j=n7+1,...n-k-1 are all zero
      !  at xi.
285   j = l5
      do i = 1, kk1
        j = j + 1
        piv = h(i)
        if (piv == 0.) cycle
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, a1(j, 1), cos, sin)
        !  transformations to right hand side.
        call fprota(cos, sin, yi, z(j))
        if (i == kk1) exit
        i2 = 1
        i3 = i + 1
        !  transformations to left hand side.
        do i1 = i3, kk1
          i2 = i2 + 1
          call fprota(cos, sin, h(i1), a1(j, i2))
        end do

      end do

      !  add contribution of this row to the sum of squares of residual
      !  right hand sides.
      fp = fp + yi**2
    end do

    fpint(n) = fp0
    fpint(n - 1) = fpold
    nrdata(n) = nplus
    !  backward substitution to obtain the b-spline coefficients c(j),j=1,.n
    call fpbacp(a1, a2, z, n7, kk, c, kk1, nest)
    !  calculate from condition (**) the coefficients c(j+n7),j=1,2,...k.
    do i = 1, k
      j = i + n7
      c(j) = c(i)
    end do

    if (iopt < 0) go to 660
    !  test whether the approximation sinf(x) is an acceptable solution.
    fpms = fp - s
    if (abs(fpms) < acc) go to 660
    !  if f(p=inf) < s accept the choice of knots.
    if (fpms < 0.) go to 350
    !  if n=nmax, sinf(x) is an interpolating spline.
    if (n == nmax) go to 630
    !  increase the number of knots.
    !  if n=nest we cannot increase the number of knots because of the
    !  storage capacity limitation.
    if (n == nest) go to 620
    !  determine the number of knots nplus we are going to add.
    npl1 = nplus * 2
    rn = nplus
    if (fpold - fp > acc) npl1 = int(rn * fpms / (fpold - fp))
    nplus = min0(nplus * 2, max0(npl1, nplus / 2, 1))
    fpold = fp
    !  compute the sum(wi*(yi-s(xi))**2) for each knot interval
    !  t(j+k) <= xi <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
    fpart = 0.
    i = 1
    l = k1
    do it = 1, m1
      if (x(it) < t(l)) go to 300
      new = 1
      l = l + 1
300   term = 0.
      l0 = l - k2
      do j = 1, k1
        l0 = l0 + 1
        term = term + c(l0) * q(it, j)
      end do

      term = (w(it) * (term - y(it)))**2
      fpart = fpart + term
      if (new == 0) cycle
      if (l <= k2) then
        fpint(nrint) = term
        new = 0
      else
        store = term * half
        fpint(i) = fpart - store
        i = i + 1
        fpart = store
        new = 0
      end if
    end do

    fpint(nrint) = fpint(nrint) + fpart
    do l = 1, nplus
      !  add a new knot
      call fpknot(x, m, t, n, fpint, nrdata, nrint, nest, 1)
      !  if n=nmax we locate the knots as for interpolation.
      if (n == nmax) go to 5
      !  test whether we cannot further increase the number of knots.
      if (n == nest) cycle
    end do

    !  restart the computations with the new set of knots.
  end do

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  part 2: determination of the smoothing periodic spline sp(x).       c
  !  *************************************************************       c
  !  we have determined the number of knots and their position.          c
  !  we now compute the b-spline coefficients of the smoothing spline    c
  !  sp(x). the observation matrix a is extended by the rows of matrix   c
  !  b expressing that the kth derivative discontinuities of sp(x) at    c
  !  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
  !  ponding weights of these additional rows are set to 1/sqrt(p).      c
  !  iteratively we then have to determine the value of p such that      c
  !  f(p)=sum(w(i)*(y(i)-sp(x(i)))**2) be = s. we already know that      c
  !  the least-squares constant function corresponds to p=0, and that    c
  !  the least-squares periodic spline corresponds to p=infinity. the    c
  !  iteration process which is proposed here, makes use of rational     c
  !  interpolation. since f(p) is a convex and strictly decreasing       c
  !  function of p, it can be approximated by a rational function        c
  !  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
  !  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
  !  to calculate the new value of p such that r(p)=s. convergence is    c
  !  guaranteed by taking f1>0 and f3<0.                                 c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  evaluate the discontinuity jump of the kth derivative of the
  !  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
350 call fpdisc(t, n, k2, b, nest)
  !  initial value for p.
  p1 = 0.
  f1 = fp0 - s
  p3 = -one
  f3 = fpms
  n11 = n10 - 1
  n8 = n7 - 1
  p = 0.
  l = n7
  do i = 1, k
    j = k + 1 - i
    p = p + a2(l, j)
    l = l - 1
    if (l == 0) go to 356
  end do

  p = p + sum(a1(:n10, 1))
356 rn = n7
  p = rn / p
  ich1 = 0
  ich3 = 0
  !  iteration process to find the root of f(p) = s.
  do iter = 1, maxit
    !  form the matrix g  as the matrix a extended by the rows of matrix b.
    !  the rows of matrix b with weight 1/p are rotated into
    !  the triangularised observation matrix a.
    !  after triangularisation our n7 x n7 matrix g takes the form
    !            ! g1 '    !
    !        g = !    ' g2 !
    !            ! 0  '    !
    !  with g2 a n7 x (k+1) matrix and g1 a n11 x n11 upper triangular
    !  matrix of bandwidth k+2. ( n11 = n7-k-1)
    pinv = one / p
    !  store matrix a into g

    c(:n7) = z(:n7)
    g1(:n7, k1) = a1(:n7, k1)
    g1(:n7, k2) = 0._8
    g2(:n7, 1) = 0._8
    g1(:n7, :k) = a1(:n7, :k)
    g2(:n7, 2:k + 1) = a2(:n7, :k)
    l = n10
    do j = 1, k1
      if (l <= 0) go to 375
      g2(l, 1) = a1(l, j)
      l = l - 1
    end do

375 do it = 1, n8
      !  fetch a new row of matrix b and store it in the arrays h1 (the part
      !  with respect to g1) and h2 (the part with respect to g2).
      yi = 0.
      h1(:k1) = 0._8
      h1(k2) = 0._8
      h2(:k1) = 0._8

      if (it > n11) go to 420
      l = it
      l0 = it
      do j = 1, k2
        if (l0 == n10) go to 400
        h1(j) = b(it, j) * pinv
        l0 = l0 + 1
      end do

      go to 470
400   l0 = 1
      do l1 = j, k2
        h2(l0) = b(it, l1) * pinv
        l0 = l0 + 1
      end do

      go to 470
420   l = 1
      i = it - n10
      do j = 1, k2
        i = i + 1
        l0 = i
430     l1 = l0 - k1
        if (l1 <= 0) go to 450
        if (l1 <= n11) go to 440
        l0 = l1 - n11
        go to 430
440     h1(l1) = b(it, j) * pinv
        cycle
450     h2(l0) = h2(l0) + b(it, j) * pinv
      end do

      if (n11 <= 0) go to 510
      !  rotate this row into triangle by givens transformations without
      !  square roots.
      !  rotation with the rows l,l+1,...n11.
470   do j = l, n11
        piv = h1(1)
        !  calculate the parameters of the givens transformation.
        call fpgivs(piv, g1(j, 1), cos, sin)
        !  transformation to right hand side.
        call fprota(cos, sin, yi, c(j))
        !  transformation to the left hand side with respect to g2.
        do i = 1, k1
          call fprota(cos, sin, h2(i), g2(j, i))
        end do

        if (j == n11) exit
        i2 = min0(n11 - j, k1)
        !  transformation to the left hand side with respect to g1.
        do i = 1, i2
          i1 = i + 1
          call fprota(cos, sin, h1(i1), g1(j, i1))
          h1(i) = h1(i1)
        end do

        h1(i1) = 0._8
      end do

      !  rotation with the rows n11+1,...n7
510   do j = 1, k1
        ij = n11 + j
        if (ij <= 0) cycle
        piv = h2(j)
        !  calculate the parameters of the givens transformation
        call fpgivs(piv, g2(ij, j), cos, sin)
        !  transformation to the right hand side.
        call fprota(cos, sin, yi, c(ij))
        if (j == k1) exit
        j1 = j + 1
        !  transformation to the left hand side.
        do i = j1, k1
          call fprota(cos, sin, h2(i), g2(ij, i))
        end do
      end do
    end do

    !  backward substitution to obtain the b-spline coefficients
    !  c(j),j=1,2,...n7 of sp(x).
    call fpbacp(g1, g2, c, n7, k1, c, k2, nest)
    !  calculate from condition (**) the b-spline coefficients c(n7+j),j=1,.
    do i = 1, k
      j = i + n7
      c(j) = c(i)
    end do

    !  computation of f(p).
    fp = 0.
    l = k1
    do it = 1, m1
      if (x(it) >= t(l)) l = l + 1
      l0 = l - k2
      term = 0.
      do j = 1, k1
        l0 = l0 + 1
        term = term + c(l0) * q(it, j)
      end do

      fp = fp + (w(it) * (term - y(it)))**2
    end do

    !  test whether the approximation sp(x) is an acceptable solution.
    fpms = fp - s
    if (abs(fpms) < acc) go to 660
    !  test whether the maximal number of iterations is reached.
    if (iter == maxit) go to 600
    !  carry out one more step of the iteration process.
    p2 = p
    f2 = fpms
    if (ich3 /= 0) go to 580
    if ((f2 - f3) > acc) go to 575
    !  our initial choice of p is too large.
    p3 = p2
    f3 = f2
    p = p * con4
    if (p <= p1) p = p1 * con9 + p2 * con1
    cycle
575 if (f2 < 0.) ich3 = 1
580 if (ich1 /= 0) go to 590
    if ((f1 - f2) > acc) go to 585
    !  our initial choice of p is too small
    p1 = p2
    f1 = f2
    p = p / con4
    if (p3 < 0.) cycle
    if (p >= p3) p = p2 * con1 + p3 * con9
    cycle
585 if (f2 > 0.) ich1 = 1
    !  test whether the iteration process proceeds as theoretically
    !  expected.
590 if (f2 >= f1 .or. f2 <= f3) go to 610
    !  find the new value for p.
    p = fprati(p1, f1, p2, f2, p3, f3)

  end do

  !  error codes and messages.
600 ier = 3
  go to 660
610 ier = 2
  go to 660
620 ier = 1
  go to 660
630 ier = -1
  go to 660
640 ier = -2
  !  the least-squares constant function c1 is a solution of our problem.
  !  a constant function is a spline of degree k with all b-spline
  !  coefficients equal to that constant c1.
  do i = 1, k1
    rn = k1 - i
    t(i) = x(1) - rn * per
    c(i) = c1
    j = i + k1
    rn = i - 1
    t(j) = x(m) + rn * per
  end do

  n = nmin
  fp = fp0
  fpint(n) = fp0
  fpint(n - 1) = 0.
  nrdata(n) = 0
660 return
end subroutine fpperi

real(8) function fprati(p1, f1, p2, f2, p3, f3)
  !  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
  !  gives the value of p such that the rational interpolating function
  !  of the form r(p) = (u*p+v)/(p+w) equals zero at p.
  !  ..
  !  ..scalar arguments..
  real(8) :: p1, f1, p2, f2, p3, f3
  !  ..local scalars..
  real(8) :: h1, h2, h3, p
  !  ..
  if (p3 <= 0._8) then
    !  value of p in case p3 = infinity.
    p = (p1 * (f1 - f3) * f2 - p2 * (f2 - f3) * f1) / ((f1 - f2) * f3)
  else
    !  value of p in case p3 ^= infinity.
    h1 = f1 * (f2 - f3)
    h2 = f2 * (f3 - f1)
    h3 = f3 * (f1 - f2)
    p = -(p1 * p2 * h3 + p2 * p3 * h1 + p3 * p1 * h2) / (p1 * h1 + p2 * h2 + p3 * h3)
  end if

  !  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
  if (f2 >= 0._8) then
    p1 = p2
    f1 = f2
  else
    p3 = p2
    f3 = f2
  end if
  fprati = p
end function fprati

subroutine fprota(cos, sin, a, b)
  !  subroutine fprota applies a givens rotation to a and b.
  !  ..
  !  ..scalar arguments..
  real(8) :: cos, sin, a, b
  ! ..local scalars..
  real(8) :: stor1, stor2
  !  ..
  stor1 = a
  stor2 = b
  b = cos * stor2 + sin * stor1
  a = cos * stor1 - sin * stor2
end subroutine fprota

subroutine fpcuro(a, b, c, d, x, n)
  !  subroutine fpcuro finds the real zeros of a cubic polynomial
  !  p(x) = a*x**3+b*x**2+c*x+d.
  !
  !  calling sequence:
  !     call fpcuro(a,b,c,d,x,n)
  !
  !  input parameters:
  !    a,b,c,d: real values, containing the coefficients of p(x).
  !
  !  output parameters:
  !    x      : real array,length 3, which contains the real zeros of p(x)
  !    n      : integer, giving the number of real zeros of p(x).
  !  ..
  !  ..scalar arguments..
  real(8) :: a, b, c, d
  integer :: n
  !  ..array argument..
  real(8) :: x(3)
  !  ..local scalars..
  integer :: i
  real(8) :: a1, b1, c1, df, disc, d1, e3, f, four, half, ovfl, pi3, p3, q, r
  real(8) :: step, tent, three, two, u, u1, u2, y
  !  ..function references..
  ! real(8) :: abs, max, datan, atan2, cos, sign, sqrt
  !  set constants
  two = 2._8
  three = 3._8
  four = 4._8
  ovfl = 1.e+04_8
  half = 0.5_8
  tent = 0.1_8
  e3 = tent / 0.3_8
  pi3 = atan(1._8) / 0.75_8
  a1 = abs(a)
  b1 = abs(b)
  c1 = abs(c)
  d1 = abs(d)

  if (max(b1, c1, d1) < a1 * ovfl) then
    !  p(x) is a third degree polynomial.
    b1 = b / a * e3
    c1 = c / a
    d1 = d / a
    q = c1 * e3 - b1 * b1
    r = b1 * b1 * b1 + (d1 - b1 * c1) * half
    disc = q * q * q + r * r
    if (disc > 0._8) then
      u = sqrt(disc)
      u1 = -r + u
      u2 = -r - u
      n = 1
      x(1) = sign(abs(u1)**e3, u1) + sign(abs(u2)**e3, u2) - b1
    else
      u = sqrt(abs(q))
      if (r < 0.) u = -u
      p3 = atan2(sqrt(-disc), abs(r)) * e3
      u2 = u + u
      n = 3
      x(1) = -u2 * cos(p3) - b1
      x(2) = u2 * cos(pi3 - p3) - b1
      x(3) = u2 * cos(pi3 + p3) - b1
    end if

  else if (max(c1, d1) < b1 * ovfl) then
    !   p(x) is a second degree polynomial.
    disc = c * c - four * b * d
    n = 0
    if (disc < 0._8) return
    n = 2
    u = sqrt(disc)
    b1 = b + b
    x(1) = (-c + u) / b1
    x(2) = (-c - u) / b1

  else if (d1 < c1 * ovfl) then
    !   p(x) is a first degree polynomial.
    n = 1
    x(1) = -d / c
  else
    !  p(x) is a constant function.
    n = 0
    return
  end if

  !
  !  apply a newton iteration to improve the accuracy of the roots.
  do i = 1, n
    y = x(i)
    f = ((a * y + b) * y + c) * y + d
    df = (three * a * y + two * b) * y + c
    step = 0.
    if (abs(f) < abs(df) * tent) step = f / df
    x(i) = y - step
  end do
end subroutine fpcuro

subroutine fpintb(t, n, bint, nk1, x, y)
  implicit none
  !  subroutine fpintb calculates integrals of the normalized b-splines
  !  nj,k+1(x) of degree k, defined on the set of knots t(j),j=1,2,...n.
  !  it makes use of the formulae of gaffney for the calculation of
  !  indefinite integrals of b-splines.
  !
  !  calling sequence:
  !     call fpintb(t,n,bint,nk1,x,y)
  !
  !  input parameters:
  !    t    : real array,length n, containing the position of the knots.
  !    n    : integer value, giving the number of knots.
  !    nk1  : integer value, giving the number of b-splines of degree k,
  !           defined on the set of knots ,i.e. nk1 = n-k-1.
  !    x,y  : real values, containing the end points of the integration
  !           interval.
  !  output parameter:
  !    bint : array,length nk1, containing the integrals of the b-splines.
  !  ..
  !  ..scalars arguments..
  integer :: n, nk1
  real(8) :: x, y
  !  ..array arguments..
  real(8) :: t(n), bint(nk1)
  !  ..local scalars..
  integer :: i, ia, ib, it, j, j1, k, k1, l, li, lj, lk
  real(8) :: a, ak, arg, b, f, one
  !  ..local arrays..
  real(8) :: aint(6), h(6), h1(6)
  logical :: minus
  !  initialization.
  one = 1.-8
  k1 = n - nk1
  ak = k1
  k = k1 - 1

  bint(:nk1) = 0._8

  !  the integration limits are arranged in increasing order.
  if (x < y) then
    a = x; b = y
    minus = .False.
  else if (x > y) then
    a = y; b = x
    minus = .True.
  else
    return
  end if

  if (a < t(k1)) a = t(k1)
  if (b > t(nk1 + 1)) b = t(nk1 + 1)
  if (a > b) return
  !  using the expression of gaffney for the indefinite integral of a
  !  b-spline we find that
  !  bint(j) = (t(j+k+1)-t(j))*(res(j,b)-res(j,a))/(k+1)
  !    where for t(l) <= x < t(l+1)
  !    res(j,x) = 0, j=1,2,...,l-k-1
  !             = 1, j=l+1,l+2,...,nk1
  !             = aint(j+k-l+1), j=l-k,l-k+1,...,l
  !               = sumi((x-t(j+i))*nj+i,k+1-i(x)/(t(j+k+1)-t(j+i)))
  !                 i=0,1,...,k
  l = k1

  !  set arg = a.
  arg = a
  do it = 1, 2
    !  search for the knot interval t(l) <= arg < t(l+1).
    do while (arg >= t(l + 1) .and. l /= nk1)
      l = l + 1
    end do

    !  calculation of aint(j), j=1,2,...,k+1.
    !  initialization.

    aint(:k1) = 0._8
    aint(1) = (arg - t(l)) / (t(l + 1) - t(l))
    h1(1) = one
    do j = 1, k
      !  evaluation of the non-zero b-splines of degree j at arg,i.e.
      !    h(i+1) = nl-j+i,j(arg), i=0,1,...,j.
      h(1) = 0._8
      do i = 1, j
        li = l + i
        lj = li - j
        f = h1(i) / (t(li) - t(lj))
        h(i) = h(i) + f * (t(li) - arg)
        h(i + 1) = f * (arg - t(lj))
      end do

      !  updating of the integrals aint.
      j1 = j + 1
      do i = 1, j1
        li = l + i
        lj = li - j1
        aint(i) = aint(i) + h(i) * (arg - t(lj)) / (t(li) - t(lj))
        h1(i) = h(i)
      end do

    end do

    if (it /= 2) exit
    !  updating of the integrals bint
    lk = l - k
    ia = lk
    do i = 1, k1
      bint(lk) = -aint(i)
      lk = lk + 1
    end do
    arg = b
  end do

  !  updating of the integrals bint.
  lk = l - k
  ib = lk - 1
  do i = 1, k1
    bint(lk) = bint(lk) + aint(i)
    lk = lk + 1
  end do

  if (ib >= ia) bint(ia:ib) = bint(ia:ib) + one

  !  the scaling factors are taken into account.
  f = one / ak
  do i = 1, nk1
    j = i + k1
    bint(i) = bint(i) * (t(j) - t(i)) * f
  end do

  !  the order of the integration limits is taken into account.
  IF (minus) bint(:nk1) = -bint(:nk1)

end subroutine fpintb

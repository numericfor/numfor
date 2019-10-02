subroutine parcur(iopt, ipar, idim, m, u, mx, x, w, ub, ue, k, s, nest, n, t, nc, c, fp, wrk, lwrk, iwrk, ier)
  !  given the ordered set of m points x(i) in the idim-dimensional space
  !  and given also a corresponding set of strictly increasing values u(i)
  !  and the set of positive numbers w(i),i=1,2,...,m, subroutine parcur
  !  determines a smooth approximating spline curve s(u), i.e.
  !      x1 = s1(u)
  !      x2 = s2(u)       ub <= u <= ue
  !      .........
  !      xidim = sidim(u)
  !  with sj(u),j=1,2,...,idim spline functions of degree k with common
  !  knots t(j),j=1,2,...,n.
  !  if ipar=1 the values ub,ue and u(i),i=1,2,...,m must be supplied by
  !  the user. if ipar=0 these values are chosen automatically by parcur
  !  as  v(1) = 0
  !      v(i) = v(i-1) + dist(x(i),x(i-1)) ,i=2,3,...,m
  !      u(i) = v(i)/v(m) ,i=1,2,...,m
  !      ub = u(1) = 0, ue = u(m) = 1.
  !  if iopt=-1 parcur calculates the weighted least-squares spline curve
  !  according to a given set of knots.
  !  if iopt>=0 the number of knots of the splines sj(u) and the position
  !  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
  !  ness of s(u) is then achieved by minimalizing the discontinuity
  !  jumps of the k-th derivative of s(u) at the knots t(j),j=k+2,k+3,...,
  !  n-k-1. the amount of smoothness is determined by the condition that
  !  f(p)=sum((w(i)*dist(x(i),s(u(i))))**2) be <= s, with s a given non-
  !  negative constant, called the smoothing factor.
  !  the fit s(u) is given in the b-spline representation and can be
  !  evaluated by means of subroutine curev.
  !
  !  calling sequence:
  !     call parcur(iopt,ipar,idim,m,u,mx,x,w,ub,ue,k,s,nest,n,t,nc,c,
  !    * fp,wrk,lwrk,iwrk,ier)
  !
  !  parameters:
  !   iopt  : integer flag. on entry iopt must specify whether a weighted
  !           least-squares spline curve (iopt=-1) or a smoothing spline
  !           curve (iopt=0 or 1) must be determined.if iopt=0 the routine
  !           will start with an initial set of knots t(i)=ub,t(i+k+1)=ue,
  !           i=1,2,...,k+1. if iopt=1 the routine will continue with the
  !           knots found at the last call of the routine.
  !           attention: a call with iopt=1 must always be immediately
  !           preceded by another call with iopt=1 or iopt=0.
  !           unchanged on exit.
  !   ipar  : integer flag. on entry ipar must specify whether (ipar=1)
  !           the user will supply the parameter values u(i),ub and ue
  !           or whether (ipar=0) these values are to be calculated by
  !           parcur. unchanged on exit.
  !   idim  : integer. on entry idim must specify the dimension of the
  !           curve. 0 < idim < 11.
  !           unchanged on exit.
  !   m     : integer. on entry m must specify the number of data points.
  !           m > k. unchanged on exit.
  !   u     : real array of dimension at least (m). in case ipar=1,before
  !           entry, u(i) must be set to the i-th value of the parameter
  !           variable u for i=1,2,...,m. these values must then be
  !           supplied in strictly ascending order and will be unchanged
  !           on exit. in case ipar=0, on exit,array u will contain the
  !           values u(i) as determined by parcur.
  !   mx    : integer. on entry mx must specify the actual dimension of
  !           the array x as declared in the calling (sub)program. mx must
  !           not be too small (see x). unchanged on exit.
  !   x     : real array of dimension at least idim*m.
  !           before entry, x(idim*(i-1)+j) must contain the j-th coord-
  !           inate of the i-th data point for i=1,2,...,m and j=1,2,...,
  !           idim. unchanged on exit.
  !   w     : real array of dimension at least (m). before entry, w(i)
  !           must be set to the i-th value in the set of weights. the
  !           w(i) must be strictly positive. unchanged on exit.
  !           see also further comments.
  !   ub,ue : real values. on entry (in case ipar=1) ub and ue must
  !           contain the lower and upper bound for the parameter u.
  !           ub <=u(1), ue>= u(m). if ipar = 0 these values will
  !           automatically be set to 0 and 1 by parcur.
  !   k     : integer. on entry k must specify the degree of the splines.
  !           1<=k<=5. it is recommended to use cubic splines (k=3).
  !           the user is strongly dissuaded from choosing k even,together
  !           with a small s-value. unchanged on exit.
  !   s     : real.on entry (in case iopt>=0) s must specify the smoothing
  !           factor. s >=0. unchanged on exit.
  !           for advice on the choice of s see further comments.
  !   nest  : integer. on entry nest must contain an over-estimate of the
  !           total number of knots of the splines returned, to indicate
  !           the storage space available to the routine. nest >=2*k+2.
  !           in most practical situation nest=m/2 will be sufficient.
  !           always large enough is nest=m+k+1, the number of knots
  !           needed for interpolation (s=0). unchanged on exit.
  !   n     : integer.
  !           unless ier = 10 (in case iopt >=0), n will contain the
  !           total number of knots of the smoothing spline curve returned
  !           if the computation mode iopt=1 is used this value of n
  !           should be left unchanged between subsequent calls.
  !           in case iopt=-1, the value of n must be specified on entry.
  !   t     : real array of dimension at least (nest).
  !           on successful exit, this array will contain the knots of the
  !           spline curve,i.e. the position of the interior knots t(k+2),
  !           t(k+3),..,t(n-k-1) as well as the position of the additional
  !           t(1)=t(2)=...=t(k+1)=ub and t(n-k)=...=t(n)=ue needed for
  !           the b-spline representation.
  !           if the computation mode iopt=1 is used, the values of t(1),
  !           t(2),...,t(n) should be left unchanged between subsequent
  !           calls. if the computation mode iopt=-1 is used, the values
  !           t(k+2),...,t(n-k-1) must be supplied by the user, before
  !           entry. see also the restrictions (ier=10).
  !   nc    : integer. on entry nc must specify the actual dimension of
  !           the array c as declared in the calling (sub)program. nc
  !           must not be too small (see c). unchanged on exit.
  !   c     : real array of dimension at least (nest*idim).
  !           on successful exit, this array will contain the coefficients
  !           in the b-spline representation of the spline curve s(u),i.e.
  !           the b-spline coefficients of the spline sj(u) will be given
  !           in c(n*(j-1)+i),i=1,2,...,n-k-1 for j=1,2,...,idim.
  !   fp    : real. unless ier = 10, fp contains the weighted sum of
  !           squared residuals of the spline curve returned.
  !   wrk   : real array of dimension at least m*(k+1)+nest*(6+idim+3*k).
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
  !    ier=0  : normal return. the curve returned has a residual sum of
  !             squares fp such that abs(fp-s)/s <= tol with tol a relat-
  !             ive tolerance set to 0.001 by the program.
  !    ier=-1 : normal return. the curve returned is an interpolating
  !             spline curve (fp=0).
  !    ier=-2 : normal return. the curve returned is the weighted least-
  !             squares polynomial curve of degree k.in this extreme case
  !             fp gives the upper bound fp0 for the smoothing factor s.
  !    ier=1  : error. the required storage space exceeds the available
  !             storage space, as specified by the parameter nest.
  !             probably causes : nest too small. if nest is already
  !             large (say nest > m/2), it may also indicate that s is
  !             too small
  !             the approximation returned is the least-squares spline
  !             curve according to the knots t(1),t(2),...,t(n). (n=nest)
  !             the parameter fp gives the corresponding weighted sum of
  !             squared residuals (fp>s).
  !    ier=2  : error. a theoretically impossible result was found during
  !             the iteration process for finding a smoothing spline curve
  !             with fp = s. probably causes : s too small.
  !             there is an approximation returned but the corresponding
  !             weighted sum of squared residuals does not satisfy the
  !             condition abs(fp-s)/s < tol.
  !    ier=3  : error. the maximal number of iterations maxit (set to 20
  !             by the program) allowed for finding a smoothing curve
  !             with fp=s has been reached. probably causes : s too small
  !             there is an approximation returned but the corresponding
  !             weighted sum of squared residuals does not satisfy the
  !             condition abs(fp-s)/s < tol.
  !    ier=10 : error. on entry, the input data are controlled on validity
  !             the following restrictions must be satisfied.
  !             -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
  !             0<=ipar<=1, 0<idim<=10, lwrk>=(k+1)*m+nest*(6+idim+3*k),
  !             nc>=nest*idim
  !             if ipar=0: sum j=1,idim (x(idim*i+j)-x(idim*(i-1)+j))**2>0
  !                        i=1,2,...,m-1.
  !             if ipar=1: ub<=u(1)<u(2)<...<u(m)<=ue
  !             if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
  !                         ub<t(k+2)<t(k+3)<...<t(n-k-1)<ue
  !                            (ub=0 and ue=1 in case ipar=0)
  !                       the schoenberg-whitney conditions, i.e. there
  !                       must be a subset of data points uu(j) such that
  !                         t(j) < uu(j) < t(j+k+1), j=1,2,...,n-k-1
  !             if iopt>=0: s>=0
  !                         if s=0 : nest >= m+k+1
  !             if one of these conditions is found to be violated,control
  !             is immediately repassed to the calling program. in that
  !             case there is no approximation returned.
  !
  !  further comments:
  !   by means of the parameter s, the user can control the tradeoff
  !   between closeness of fit and smoothness of fit of the approximation.
  !   if s is too large, the curve will be too smooth and signal will be
  !   lost ; if s is too small the curve will pick up too much noise. in
  !   the extreme cases the program will return an interpolating curve if
  !   s=0 and the least-squares polynomial curve of degree k if s is
  !   very large. between these extremes, a properly chosen s will result
  !   in a good compromise between closeness of fit and smoothness of fit.
  !   to decide whether an approximation, corresponding to a certain s is
  !   satisfactory the user is highly recommended to inspect the fits
  !   graphically.
  !   recommended values for s depend on the weights w(i). if these are
  !   taken as 1/d(i) with d(i) an estimate of the standard deviation of
  !   x(i), a good s-value should be found in the range (m-sqrt(2*m),m+
  !   sqrt(2*m)). if nothing is known about the statistical error in x(i)
  !   each w(i) can be set equal to one and s determined by trial and
  !   error, taking account of the comments above. the best is then to
  !   start with a very large value of s ( to determine the least-squares
  !   polynomial curve and the upper bound fp0 for s) and then to
  !   progressively decrease the value of s ( say by a factor 10 in the
  !   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
  !   approximating curve shows more detail) to obtain closer fits.
  !   to economize the search for a good s-value the program provides with
  !   different modes of computation. at the first call of the routine, or
  !   whenever he wants to restart with the initial set of knots the user
  !   must set iopt=0.
  !   if iopt=1 the program will continue with the set of knots found at
  !   the last call of the routine. this will save a lot of computation
  !   time if parcur is called repeatedly for different values of s.
  !   the number of knots of the spline returned and their location will
  !   depend on the value of s and on the complexity of the shape of the
  !   curve underlying the data. but, if the computation mode iopt=1 is
  !   used, the knots returned may also depend on the s-values at previous
  !   calls (if these were smaller). therefore, if after a number of
  !   trials with different s-values and iopt=1, the user can finally
  !   accept a fit as satisfactory, it may be worthwhile for him to call
  !   parcur once more with the selected value for s but now with iopt=0.
  !   indeed, parcur may then return an approximation of the same quality
  !   of fit but with fewer knots and therefore better if data reduction
  !   is also an important objective for the user.
  !
  !   the form of the approximating curve can strongly be affected by
  !   the choice of the parameter values u(i). if there is no physical
  !   reason for choosing a particular parameter u, often good results
  !   will be obtained with the choice of parcur (in case ipar=0), i.e.
  !        v(1)=0, v(i)=v(i-1)+q(i), i=2,...,m, u(i)=v(i)/v(m), i=1,..,m
  !   where
  !        q(i)= sqrt(sum j=1,idim (xj(i)-xj(i-1))**2 )
  !   other possibilities for q(i) are
  !        q(i)= sum j=1,idim (xj(i)-xj(i-1))**2
  !        q(i)= sum j=1,idim abs(xj(i)-xj(i-1))
  !        q(i)= max j=1,idim abs(xj(i)-xj(i-1))
  !        q(i)= 1
  !
  !  other subroutines required:
  !    fpback,fpbspl,fpchec,fppara,fpdisc,fpgivs,fpknot,fprati,fprota
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
  implicit none
  integer, intent(IN) :: iopt !< Flag. Possible values (-1, 0, 1)
  integer, intent(IN) :: ipar !< Flag. Possible values (0, 1)
  integer, intent(IN) :: idim !< Dimension of curve (0 < idim < 11)
  integer, intent(IN) :: m !< size(u)
  integer, intent(IN) :: mx !< size(x)
  integer, intent(IN) :: k !< Degree of spline
  integer, intent(IN) :: nest !< Must be >= that number of knots
  integer, intent(INOUT) :: n !< Will have the total number of knots
  integer, intent(IN) :: nc !< size(c)
  integer, intent(OUT) :: ier !< Flag. Status of output
  integer, intent(IN) :: lwrk !< size(wrk)
  !  ..array arguments..
  integer, dimension(nest), intent(INOUT) :: iwrk !< integer workspace
  real(8), intent(INOUT) :: ub !< lower bound
  real(8), intent(INOUT) :: ue !< upper bound
  real(8), intent(IN) :: s !< Smoothing parameter
  real(8), intent(OUT) :: fp !< squared residuals
  real(8), dimension(m), intent(INOUT) :: u !<
  real(8), dimension(mx), intent(IN) :: x !< Array of size = idim*m
  real(8), dimension(m), intent(IN) :: w !< Weights (size >= m)
  real(8), dimension(nest), intent(INOUT) :: t !< Knots.
  real(8), dimension(nc), intent(OUT) :: c !< Array of size >= nest*idim.
  real(8), dimension(lwrk), intent(INOUT) :: wrk !<

  !  ..local scalars..
  real(8) :: tol, dist
  integer :: i, ia, ib, ifp, ig, iq, iz, i1, i2, j, k1, k2, lwest, maxit, nmin, ncc
  ! ..function references
  ! real(8) :: sqrt
  !  ..
  !  we set up the parameters tol and maxit
  maxit = 20
  tol = 0.001_8
  !  before starting computations data checks are made. If the input data
  !  are invalid, control is immediately repassed to the calling program.
  ier = 10
  IF (iopt < (-1) .or. iopt > 1) return ! Check the task flag
  IF (ipar < 0 .or. ipar > 1) return ! Check the boundary values flag
  IF (idim <= 0 .or. idim > 10) return ! Check correct dimension of curve
  IF (k <= 0 .or. k > 5) return ! Check the order of the splines
  IF (any(w <= 0._8)) return    ! Check the weigths

  k1 = k + 1
  k2 = k1 + 1
  nmin = 2 * k1
  IF (m < k1 .or. nest < nmin) return
  ncc = nest * idim
  IF (mx < m * idim .or. nc < ncc) return
  lwest = m * k1 + nest * (6 + idim + 3 * k)
  IF (lwrk < lwest) return

  if (ipar == 0 .and. iopt <= 0) then ! Determine u(1:m) automatically
    u(1) = 0._8
    do i = 2, m
      dist = sqrt(sum((x(idim * (i - 1) + 1:idim * i) - x(idim * (i - 2) + 1:idim * (i - 1)))**2))
      u(i) = u(i - 1) + dist
    end do

    IF (u(m) <= 0._8) return
    u(2:m) = u(2:m) / u(m)

    ub = 0._8
    ue = 1._8
    u(m) = ue
  end if

  IF (ub > u(1) .or. ue < u(m)) return
  IF (any(u(1:m - 1) >= u(2:m))) return

  if (iopt == -1) then
    IF (n < nmin .or. n > nest) return
    t(1:k1) = ub
    t(n - k:n) = ue
    call fpchec(u, m, t, n, k, ier)
    IF (ier /= 0) return
  else
    IF (s < 0._8) return ! JF: This should be at the beginning of the routine?
    IF (s == 0._8 .and. nest < (m + k1)) return
    ier = 0
  end if

  ! We partition the working space and determine the spline curve.
  ifp = 1
  iz = ifp + nest
  ia = iz + ncc
  ib = ia + nest * k1
  ig = ib + nest * k2
  iq = ig + nest * k2
  call fppara(iopt, idim, m, u, mx, x, w, ub, ue, k, s, nest, tol, maxit, k1, k2,&
    &n, t, ncc, c, fp, wrk(ifp), wrk(iz), wrk(ia), wrk(ib), wrk(ig), wrk(iq),&
    &iwrk, ier)
end subroutine parcur

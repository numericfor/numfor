!! @date "2019-11-21 11:49:03"

!> quads is a wrapper for simple use of a subset of quadpack integrators
!! with no weight function
!!
!! Examples:
!!
   ! PRNM(quad) (f, a, b, IntVal, abserr, epsabs, epsrel, maxsub, points, neval, ier, info)
   ! implicit none
   ! NUMFOR_KINDF :: f !<
   ! NUMFOR_KINDR, intent(OUT) :: IntVal !<
   ! TYPE(QPACK_EXTRA), optional, intent(INOUT) :: info !< Packs extra info from the integration routines. Default size of arrays is 500

   ! real(dp), intent(IN) :: a !<
   ! real(dp), intent(IN) :: b !<
   ! real(dp), intent(OUT) :: abserr !<
   ! real(dp), optional, intent(IN) :: epsabs !<
   ! real(dp), optional, intent(IN) :: epsrel !<
   ! integer, optional, intent(IN) :: maxsub !<
   ! real(dp), dimension(:), optional, intent(IN) :: points !<
   ! integer, optional, intent(OUT) :: ier !<
   ! integer, optional, intent(OUT) :: neval !<

   ! real(dp) :: epsabs_, epsrel_
   ! integer :: maxsub_, ier_, neval_
   ! character(len=4) :: qkrule = 'qk21'

   epsabs_ = 1.5e-8_dp; IF (Present(epsabs)) epsabs_ = epsabs
   epsrel_ = 1.5e-7_dp; IF (Present(epsrel)) epsrel_ = epsrel
   maxsub_ = 500; IF (Present(maxsub)) maxsub_ = maxsub ! "limit" in the original routines

   ! If not initialized, use default size of maxsub
   ! IF (Present(info) .and. (info%last == 0)) info%init(500)
   IF (Present(info) .and. (info%size == 0)) info = 500

   ! If there are any break points, we use them
   if (Present(points)) then
     ! Incompatible parameters
     IF (is_inf(a) /= 0 .or. is_inf(b) /= 0) call print_msg("Infinite integration limits cannot be used with break points.")
     IF (Present(info)) call print_msg("Extra parameters not implemented with break points. Continuing")

     call qagp(f, a, b, points, IntVal, epsabs_, epsrel_, abserr, neval_, ier_)
     ! for finite limits call qags
   else if (is_inf(a) == 0 .and. is_inf(b) == 0) then      ! finite limits
     if (Present(info)) then
       call qags(f, a, b, IntVal, epsabs_, epsrel_, qkrule, abserr, neval_, ier_, info)
       info%msg = "Routine used: qags"
     else
       call qags(f, a, b, IntVal, epsabs_, epsrel_, qkrule, abserr, neval_, ier_)
     end if
   else                          ! infinite integration domain
     if (Present(info)) then
       call qagi(f, a, b, IntVal, epsabs_, epsrel_, qkrule, abserr, neval_, ier_, info)
       info%msg = "Routine used: qagi"
     else
       call qagi(f, a, b, IntVal, epsabs_, epsrel_, qkrule, abserr, neval_, ier_, info)
     end if
   end if

   IF (Present(ier)) ier = ier_
   IF (Present(neval)) neval = neval_
   ! Not implemented yet

   ! end PRNM(quad)

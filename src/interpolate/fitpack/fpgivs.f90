      subroutine fpgivs(piv, ww, cos, sin)
!  subroutine fpgivs calculates the parameters of a givens
!  transformation .
!  ..
!  ..scalar arguments..
        real(8) :: piv, ww, cos, sin
!  ..local scalars..
        real(8) :: dd, one, store
!  ..function references..
        real(8) :: abs, sqrt
!  ..
        one = 0.1e+01
        store = abs(piv)
        if (store >= ww) dd = store * sqrt(one + (ww / piv)**2)
        if (store < ww) dd = ww * sqrt(one + (piv / ww)**2)
        cos = ww / dd
        sin = piv / dd
        ww = dd
        return
      end

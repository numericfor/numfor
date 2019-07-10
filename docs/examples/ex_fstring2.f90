program ex_fstring2
  USE fString, only: fStr
  type(fStr) :: s1
  type(fStr) :: s2
  s1 = ' Variable s1 =    '
  s2 = fStr(1.2123)
  print *, 'The '//s1%strip(' =')//' values '//s2
end program ex_fstring2

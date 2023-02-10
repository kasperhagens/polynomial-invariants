## Author: Kasper Hagens <kasperhagens@kaspers-air.home>
## Created: 2022-12-23

%Substitution CHI(y, z', z) = (y-1, z'+y, z+y)
%or equivalently
%CHI(Y(1), Y(2), Y(3)) = (Y(1)-1, Y(2)+Y(1), Y(3)+Y(1))

function CHI = substitutionSumtwothree (Y)
  CHI(1) = Y(1)-1;
  CHI(2) = Y(2)+Y(1);
  CHI(3)=Y(3)+Y(1);
endfunction

## Author: Kasper Hagens
## Created: 2023-01-21
%----------------------------------------------------------------------------
%INSTALLATION SYMBOLIC PACKAGE
%Installation of the symbolic package is required.
%Download symbolic-3.0.1.tar.gz on
%https://octave.sourceforge.io/symbolic/
%and store it in the current working directory.
%Now execute the following command in the command-window of GNU Octave:
%pkg install -forge symbolic-3.0.1.tar.gz
%
%After installation, it is still needed to load the package into GNU Octave.
%This is done by executing the following command:
%pkg load symbolic
%Now the function is ready to use.
%----------------------------------------------------------------------------
%FUNCTION DESCRIPTION
%This function will compute all polynomial invariants of degree d for a
%divergence pattern, consisting of divergence vectors of length n,
%initialized by some vector symbolic vector init. The inputs required by this function are
%n: the number of dicergence variables
%d: the degree of the polynomial invariant
%chi: the divergence substitution
%init: the initialization vector, which should be of length n
%
%The output will name divergence variables by Y_11,...,Y_1n.

%INITIALIZATION VECTOR
%If we, for example, have an initialization vector (x^2-1, x*y) then we can declare it by
%syms x y;
%init = sym([x^2-1,x*y]);
%
%For a given k, it is also possible to declare a symbolic vector containing
%variables x_11, x_12, ..., x_1k. For this, use the following command
%x = sym('x',[1 k]);
%
%
%!!!IMPORTANT!!! Do not declare any syms of the shape Y_ij
%since they are already reserved for the divergence variables.

%HOW TO CALL THIS FUNCTION
%The substitution is a user-defined function which should be stored in the same
%folder as this function. So this function takes another function as argument.
%This requires to put a @ in front of the chi-argument.
%So the correct way of calling this function is
%
%invariants(n,d,@chi, init)
function Invariants = invariants (n, d, chi, init)
  %We have divergence variables Y_11, Y12 ,..., Y_1n
  Y = sym('Y',[1 n]);
  J = zeros(1,n);
  M = sym(zeros(n,1));

  % This for-loop computes, for every 0<=p<=d, the set of exponents for the monomials of degree p
  % {(j_1, ..., j_n) | sum_{i=1}^n j_i = p} and puts them together into the matrix J.
  % Simultaneously, for every p, we compute the monomials of degree p:
  % {x(1)^j_1*...*x(n)^j_n | sum_{i=1}^n j_i = p} and put them together into a matrix M.
  count = 1;
  for p=[0:d]
    %This part computes the cartesian product {0,..p}^n
    [tempctprod{1:n}] = ndgrid((repmat({0:p},1,n)){:});
    ctprod = reshape(cat(n,tempctprod{:}),[],n);
    clear tempctprod;

    l=size(ctprod,1);
    for i = [1:l]
    %Here we check, for every (j_1, ..., j_n), the condition sum_{i=1}^n j_i = p
    %and in case it is true then we store the vector in matrix J and we compute
    %Y(1)^j_1*...*Y(n)^j_n
      if ( sum(ctprod(i, :)) == p )
        J(count,:) = ctprod(i, :);
        M(count)= prod(Y.^ctprod(i, :));
        count = count + 1;
      endif
    endfor
  endfor

  %The first vector of the divergence space is given by init. This vector is stored
  %in a matrix div.
  div = sym([1:n]);
  div(1, :) = init;

  %The number of monomials in the polynomial invariant equals rows(M).
  %Therefore we need to compute the first rows(M) vectors in divergence space.
  %Each succeeding vector is obtained by applying the substitution chi to the previous one.
  for i=[2:rows(M)]
        div(i, :) = chi(div(i-1, :));
  end

  %We need to use these div(i, :), for 1<=i<=n+1 as inputs for the function
  %and store them in the divergence matrix D.
  D = sym(zeros(1,rows(M)));

  for i=[1:rows(M)]
    D(i,:) = subs(M,Y,div(i, :))';
  endfor

%-----------------------------------------------------------------------------
  strcat(puts("We compute the polynomial invariants of degree "),
  puts(num2str(d)),
  puts(" corresponding to the divergence pattern initialized by "),
  disp(init),
  puts("and substitution "),
  disp(chi)
  )

  disp("The monomial vector is given by \n")
  disp(M)
  disp("\n")

  disp("The divergence matrix is given by \n")
  disp(D)
  disp("\n")

  %Compute a basis for ker(D)
  disp("A basis for its kernel is given by \n")
  B=null(D);
  disp(B)
  disp("\n")

  %Convert to invariants by taking the innerproduct with the monomial vector M.
  Inv=sym([]);
  Inv = B'*M;

  disp("corresponding to the following invariants \n")
  Invariants = sym((Inv == sym(zeros(columns(B),1))))

endfunction

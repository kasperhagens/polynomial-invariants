function Invariants = invariants (n, d, chi, init)
  %Generate divergence variables Y_11, Y12 ,..., Y_1n
  Y = sym('Y',[1 n]);
  J = zeros(1,n);
  M = sym(zeros(n,1));

  % Generate for every 0<=p<=d, the set of exponents for the monomials of degree p: 
  % {(j_1, ..., j_n) | sum_{i=1}^n j_i = p} and store them in the matrix J.
  % Simultaneously, for every p, compute the monomials of degree p:
  % {x(1)^j_1*...*x(n)^j_n | sum_{i=1}^n j_i = p} and store them in a matrix M.
  count = 1;
  for p=[0:d]
    %Generate the cartesian product {0,..p}^n
    [tempctprod{1:n}] = ndgrid((repmat({0:p},1,n)){:});
    ctprod = reshape(cat(n,tempctprod{:}),[],n);
    clear tempctprod;

    l=size(ctprod,1);
    for i = [1:l]
    % Check which exponent vectors (j_1, ..., j_n) satisfy
    % sum_{i=1}^n j_i = p, store these vectors in J, and compute the
    % corresponding monomials Y_1^j_1 * ... * Y_n^j_n.
      if ( sum(ctprod(i, :)) == p )
        J(count,:) = ctprod(i, :);
        M(count)= prod(Y.^ctprod(i, :));
        count = count + 1;
      endif
    endfor
  endfor

  %The first vector of the divergence space is given by init. This vector is stored in a matrix div.
  div = sym([1:n]);
  div(1, :) = init;

  %The number of monomials in the polynomial invariant equals rows(M).
  %Therefore, according to our heuristics, we compute the first rows(M) vectors in divergence space.
  %Each succeeding vector is obtained by applying the substitution chi to the previous one.
  for i=[2:rows(M)]
        div(i, :) = chi(div(i-1, :));
  end

  % Construct the divergence matrix D by evaluating the monomial vector M at each of the generated divergence vectors.
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

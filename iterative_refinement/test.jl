using LinearAlgebra

function ss()
# n= 1
# A = [10000000*randn(n,n) randn(n,n);
#      randn(n,n) 0.0000001*randn(n,n)]
# A = [3.9 1.6;6.8 2.9]
# A[2,:] = (ones(2) + 1e-10*randn(2)) .* A[1,:]
A = float([64       -2016       20160      -92400      221760     -288288      192192      -51480;
-2016       84672     -952560     4656960   -11642400    15567552   -10594584     2882880 ;
20160     -952560    11430720   -58212000   149688000  -204324120   141261120   -38918880 ;
-92400     4656960   -58212000   304920000  -800415000  1109908800  -776936160   216216000 ;
221760   -11642400   149688000  -800415000  2134440000 -2996753760  2118916800  -594594000 ;
-288288    15567552  -204324120  1109908800 -2996753760  4249941696 -3030051024   856215360 ;
192192   -10594584   141261120  -776936160  2118916800 -3030051024  2175421248  -618377760 ;
-51480     2882880   -38918880   216216000  -594594000   856215360  -618377760   176679360 ])

b = zeros(8)
b[3] = 1

x = A\b

@show A*x - b

@show cond(A)

println("Naive Guess")
@show norm(A*x - b)

Af = factorize(A)
# r = zeros(2*n)
for i = 1:5
    r = b-A*x
    # r .= A*x - b
    # x -= A\r
    d = Af\r
    x = x +d

    @show norm(A*x -b)
end

@show norm(A*x -b)

end

ss()

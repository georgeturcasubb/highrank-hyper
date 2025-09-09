///////////////////////////////////////////////////////////////////////
// highrankgen2.m  —  Robust script with logging to highrank.txt
// For the article:
//
//   E : Y^2 = X^3 - s1*X^2 + s2*X - 1
//   C : Y^2 = X^6 - s1*X^4 + s2*X^2 - 1
//
// Relations: u = s1*s2,  v = s1^3 + s2^3 = 9u - 27
// u(t) chosen so j(E) = j(E_t) with Kloosterman family:
//   E_t : y^2 = x^3 + (2U) x + (4 t^2 V)
//   U = t^8 + 14 t^4 + 1,   V = t^8 + 6 t^4 + 1
///////////////////////////////////////////////////////////////////////

// ---- Send everything to highrank.txt right away
SetOutputFile("highrank.txt");

// (A) Global settings
SetClassGroupBounds("GRH");  // speed-ups relying on GRH where relevant
Q := Rationals();
PRQ<x> := PolynomialRing(Q);

// Pretty banner
printf "==== highrankgen2.m started ====\n";
printf "Timestamp (Magma) : %o\n\n", Cputime();

// (B) Helpers: absolute discriminant & absolute degree to Q (robust)
function AbsDiscAndDeg(F)
    if Type(F) eq FldRat then
        return Integers()!1, 1;
    elif IsAbsoluteField(F) then
        ZF := Integers(F);
        return Abs(Discriminant(ZF)), Degree(F);
    else
        AF := AbsoluteField(F);
        ZAF := Integers(AF);
        return Abs(Discriminant(ZAF)), Degree(AF);
    end if;
end function;

procedure PrintFieldInfo(label, F)
    D, n := AbsDiscAndDeg(F);
    printf "%o: [abs deg = %o], abs discriminant = %o\n", label, n, D;
end procedure;

// (C) Robust square-root adjunction over number fields
function AdjoinSqrt(F, a)
// Returns (M, r) with M/F quadratic or equal to F, and r in M satisfying r^2 = a.
    if a eq 0 then
        return F, F!0;
    end if;
    if IsSquare(a) then                  // boolean form
        rt := SquareRoot(a);             // obtain the root explicitly
        return F, F!rt;
    else
        PF<z> := PolynomialRing(F);
        M<r> := ext<F | z^2 - a>;
        return M, r;
    end if;
end function;

// (D) Kloosterman U, V, j(E_t), and closed form u(t)
U_of_t := func< t | t^8 + 14*t^4 + 1 >;
V_of_t := func< t | t^8 + 6*t^4 + 1 >;
jEt := func< t | Q!( (3456*U_of_t(t)^3) / (2*U_of_t(t)^3 + 27*(t^4)*V_of_t(t)^2) ) >;
u_of_t := func< t | Q!(
    (9*t^24 + 864*t^20 + 11151*t^16 + 43920*t^12 + 11151*t^8 + 864*t^4 + 9) /
    (4*t^24 + 222*t^20 + 3012*t^16 + 13364*t^12 + 3012*t^8 + 222*t^4 + 4)
)>;

// (E) Choose a good t0 (or fix a specific one by replacing this block)
random_range := 80;
max_tries := 600;
good := false;  t0 := 0;  u0 := Q!0;  j0 := Q!0;

for attempt in [1..max_tries] do
    t0 := Random(-random_range, random_range);
    if t0 eq 0 then continue; end if;

    denomJ := 2*U_of_t(t0)^3 + 27*(t0^4)*V_of_t(t0)^2;
    if denomJ eq 0 then continue; end if;

    j0 := jEt(t0);
    u0 := u_of_t(t0);

    // Avoid j in {0,1728} and quadratic degeneracy 12u-27=0
    if j0 eq 0 or j0 eq 1728 then continue; end if;
    if 12*u0 - 27 eq 0 then continue; end if;

    // Guard (U(t0)!=0; never for integer t, but keep)
    if U_of_t(t0) eq 0 then continue; end if;

    good := true;  break;
end for;

if not good then
    printf "ERROR: Failed to find a good random t0; increase random_range or max_tries.\n";
    UnsetOutputFile();  quit;
end if;

printf "Chosen t0 = %o\n", t0;
printf "j(E_t)(t0) = %o\n", j0;
printf "u(t0)      = %o\n\n", u0;

// (F) Build L as splitting field of S^2 - 3S + 9 - 3u  (quadratic unless split)
Squad := x^2 - 3*x + (9 - 3*u0);
L := SplittingField(Squad);
RL<yL> := PolynomialRing(L);
rts := Roots(RL!Squad);
S := (#rts ge 1) select rts[1][1] else (L!0);

PrintFieldInfo("Field L", L);
printf "  Defining polynomial for a primitive element of L over Q:\n";
try
    printf "  %o\n\n", DefiningPolynomial(L);
catch e
    printf "  (DefiningPolynomial not available.)\n\n";
end try;

// (G) Build M = L( sqrt(S^2 - 4u) )
R2_L := L!S^2 - L!(4*u0);
M, r := AdjoinSqrt(L, R2_L);
if M eq L then
    printf "R2 = S^2 - 4u is a square in L; set M = L.\n";
else
    printf "Adjoined sqrt(S^2 - 4u) to get M.\n";
end if;

PrintFieldInfo("Field M", M);
printf "  Defining polynomial for a primitive element of M over Q:\n";
try
    printf "  %o\n\n", DefiningPolynomial(M);
catch e
    printf "  (DefiningPolynomial not available.)\n\n";
end try;

// (H) Solve s1, s2 over M
S_M := M!S;
uM  := M!u0;
s1 := (S_M + r)/2;
s2 := (S_M - r)/2;

// Sanity: u = s1*s2, v = 9u - 27
vM := s1^3 + s2^3;
assert s1*s2 eq uM;
assert vM eq M!(9*u0 - 27);

// (I) Curves over M
RM<xM> := PolynomialRing(M);
E_M := EllipticCurve([ 0, -s1, 0, s2, -1 ]);                  // y^2 = x^3 - s1 x^2 + s2 x - 1
C_M := HyperellipticCurve( xM^6 - s1*xM^4 + s2*xM^2 - 1 );    // y^2 = sextic

// (J) Build K = M( sqrt(D(t0)) ), D(t) = -(9 - u)/(6U)
U0 := U_of_t(t0);
D0 := M!( - (Q!9 - u0) / (Q!6 * U0) );
K, d := AdjoinSqrt(M, D0);
if K eq M then
    printf "D(t0) is a square in M; set K = M.\n";
else
    printf "Adjoined sqrt(D(t0)) to get K.\n";
end if;

PrintFieldInfo("Field K", K);
printf "  Defining polynomial for a primitive element of K over Q:\n";
try
    printf "  %o\n", DefiningPolynomial(K);
catch e
    printf "  (DefiningPolynomial not available.)\n";
end try;

// (K) Base-change to K
E_K := BaseChange(E_M, K);
C_K := ChangeRing(C_M, K);
J_K := Jacobian(C_K);

// RankBounds is only implemented for Jacobians over the Rationals so cannot run on J_K.

// Report j-invariants (should match j(E_t))
printf "\n--- j-invariants ---\n";
printf "j(E/K)     = %o\n", jInvariant(E_K);
printf "j(E_t)(t0) = %o\n", j0;
assert jInvariant(E_K) eq j0;

// (L) Optional Rank bounds (may be heavy; wrapped)
// printf "\n=== Rank bounds over K (with GRH class-group bounds) ===\n";
// try
//    t0 := Cputime();
//    lbE, ubE := RankBounds(E_K);
//    t1 := Cputime(t0);
//    printf "RankBounds(E_K) : [%o, %o]   [CPU time: %o seconds]\n", lbE, ubE, t1;
//catch e
//    printf "RankBounds(E_K) failed with message:\n%o\n", e`Object;
//end try;


// (M) Field summary and exit
printf "\n=== Field summary ===\n";
PrintFieldInfo("L", L);
PrintFieldInfo("M", M);
PrintFieldInfo("K", K);

printf "\n==== highrankgen2.m finished ====\n";
UnsetOutputFile();
quit;

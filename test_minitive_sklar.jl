using PossibilisticArithmetic, BivariateCopulas, IntervalArithmetic, Test

C_survival( C :: copula, u1, u2) = u1 + u2 - 1 + C(1 - u1, 1 - u2)[1]

###
#   Calculates the probabiliy measure in some 2-D rectange (box) from a 2-D cdf (joint)
###
function calc_measure( J :: Joint, box :: IntervalBox)
    lefts = left.(box)
    rights = right.(box)

    prob_measure = cdf(J, rights[1], rights[2]) - cdf(J, lefts[1], rights[2]) - cdf(J, rights[1], lefts[2]) + cdf(J, lefts[1], lefts[2])
    return prob_measure[1]
end

###
#   Calculates the belief & plausability in box from a DSS
###
function calc_bounds(focal_el :: Vector{IntervalBox{2, Float64}}, masses :: Vector{Float64}, box :: IntervalBox{2, Float64})

    belief = 0.0
    plaus  = 0.0
    for (i, el) in enumerate(focal_el)
        intersect = el ∩ box != ∅
        sub = el ⊆ box
        if intersect
            plaus += masses[i]
        end
        if sub
            belief += masses[i]
        end
    end
    return interval(belief, plaus)
end

###
#   Uses math from section 5.1 in sklar's theorem for minitive beliefs
###
function sklar_survival(x1 :: Fuzzy, x2 :: Fuzzy, C :: copula, box :: IntervalBox)
    be1 = mass(x1, box[1]).lo       # returns marginal [bel, plaus]
    be2 = mass(x2, box[2]).lo
    return C_survival(C, be1, be2)
end

###
#   Minitive sklar
###
function minitive_sklar(x1 :: Fuzzy, x2 :: Fuzzy, C :: copula, box :: IntervalBox)
    be1 = mass(x1, box[1]).lo
    be2 = mass(x2, box[2]).lo
    return C(be1, be2)[1]
end


###
# Example begins
###
a1 = Fuzzy(0, 0.5, 1, steps = 500)    # Marginal fuzzy 1
a2 = Fuzzy(0, 0.5, 1, steps = 500)    # Marginal fuzzy 2

a1Dist = Uniform(0, 0.5)               # Marginal dist 1
a2Dist = Uniform(0, 0.5)               # Marginal dist 2


@test check_inside(a1, a1Dist)      # checks if marginals are consistent
@test check_inside(a2, a2Dist)

copulas = [W(), Gaussian(-0.7), Gaussian(-0.5), Gaussian(-0.2), Pi(), Gaussian(0.2), Gaussian(0.5), Gaussian(0.7), M()]

ntests = 10^4

tolerance = 10^-14

for C in copulas

    println(" Testing: ")
    println(" F1 = $a1 ")
    println(" F2 = $a2 ")
    println(" C = $C ")

    println()
    println(" Against: ")
    println(" D1 = $a1Dist ")
    println(" D2 = $a2Dist ")
    println(" C = $C ")

    println()
    println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    println(" Only failures will be printed")
    println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    println()


    println()
    println()
    println()
    println("--------------------------------------------")
    println("checking sets:  [-∞, x] × [-∞, y]")
    println("--------------------------------------------")
    println("box                                            |      μ_xy     |      [bel, plaus]     |      bel_sklar      |      bel_survival ")


    BivDist = C(a1Dist, a2Dist)         # bivariate distribution
    masses, cart = PossibilisticArithmetic.mobiusTransform2D(a1, a2 , C)    # Should be called something else
    carts = IntervalBox.(cart)

    failures = 0

    for i = 1:ntests

        bounds = rand(2)
        box = interval(-∞, bounds[1]) × interval(-∞, bounds[2])

        prob = calc_measure(BivDist, box)
        if prob < 0; prob = 0; end

        probInt = calc_bounds(carts, masses, box)

        bel1 = sklar_survival(a1, a2, C, box)
        if bel1 < 0; bel1 = 0; end

        bel2 = minitive_sklar(a1, a2, C, box)

        t1 =  probInt.lo <= prob
        t2 = bel1 <= prob + tolerance
        t3 = bel2 <= prob
        #t4 = bel1 ≈ bel2

        if any( [!t1, !t2, !t3])
            println("$box     |     $prob     |      $(probInt)      |      $(bel2)      |      $(bel1)")
            failures += 1
        end
        #@test probInt.lo <= prob
        #@test bel1 <= prob
        #@test bel2 <= prob
        #@test bel1 ≈ bel2

    end

    println()
    println("Failures: $failures/$ntests")


    println()
    println()
    println()
    println("--------------------------------------------")
    println("checking sets:  [-∞, x] × [-∞, y]")
    println("--------------------------------------------")
    println("box                                            |      μ_xy     |      [bel, plaus]     |      bel_sklar      |      bel_survival ")

    failures = 0

    for i = 1:ntests

        bounds = rand(2)

        box = interval(-∞, bounds[1]) × interval(-∞, bounds[2])

        prob = calc_measure(BivDist, box)
        if prob < 0; prob = 0;end

        probInt = calc_bounds(carts, masses, box)

        bel1 = sklar_survival(a1, a2, C, box)
        if bel1 < 0; bel1 = 0; end
        bel2 = minitive_sklar(a1, a2, C, box)

        t1 =  probInt.lo <= prob
        t2 = bel1 <= prob + tolerance
        t3 = bel2 <= prob

        if any( [!t1, !t2, !t3])
            println("$box     |     $prob     |      $(probInt)      |      $(bel2)      |      $(bel1)")
            failures += 1
        end

    end
    println()
    println("Failures: $failures/$ntests")

    println()
    println()
    println()
    println("--------------------------------------------")
    println("checking sets:  [x, ∞] × [-∞, y]")
    println("--------------------------------------------")
    println("box                                            |      μ_xy     |      [bel, plaus]     |      bel_sklar      |      bel_survival ")

    failures = 0

    for i = 1:ntests

        bounds = rand(2)

        box = interval(bounds[1], ∞) × interval( -∞, bounds[2])

        prob = calc_measure(BivDist, box)
        if prob < 0; prob = 0;end

        probInt = calc_bounds(carts, masses, box)

        bel1 = sklar_survival(a1, a2, C, box)
        if bel1 < 0; bel1 = 0; end
        bel2 = minitive_sklar(a1, a2, C, box)

        t1 =  probInt.lo <= prob
        t2 = bel1 <= prob + tolerance
        t3 = bel2 <= prob

        if any( [!t1, !t2, !t3])
            println("$box     |     $prob     |      $(probInt)      |      $(bel2)      |      $(bel1)")
            failures += 1
        end

    end
    println()
    println("Failures: $failures/$ntests")



    println()
    println()
    println()
    println("--------------------------------------------")
    println("checking sets:  [∞, x] × [y, ∞]")
    println("--------------------------------------------")
    println("box                                            |      μ_xy     |      [bel, plaus]     |      bel_sklar      |      bel_survival ")

    failures = 0

    for i = 1:ntests

        bounds = rand(2)

        box = interval(-∞, bounds[1]) × interval(bounds[2],  ∞)

        prob = calc_measure(BivDist, box)
        if prob < 0; prob = 0;end

        probInt = calc_bounds(carts, masses, box)

        bel1 = sklar_survival(a1, a2, C, box)
        if bel1 < 0; bel1 = 0; end
        bel2 = minitive_sklar(a1, a2, C, box)

        t1 =  probInt.lo <= prob
        t2 = bel1 <= prob + tolerance
        t3 = bel2 <= prob

        if any( [!t1, !t2, !t3])
            println("$box     |     $prob     |      $(probInt)      |      $(bel2)      |      $(bel1)")
            failures += 1
        end

    end
    println()
    println("Failures: $failures/$ntests")

    println()
    println()
    println()
    println("--------------------------------------------")
    println("checking sets:  [x1, x2] × [y1, y2]")
    println("--------------------------------------------")
    println("box                                            |      μ_xy     |      [bel, plaus]     |      bel_sklar      |      bel_survival ")

    failures = 0

    for i = 1:ntests

        bounds = rand(4)

        b1 = sort(bounds[1:2])
        b2 = sort(bounds[3:4])

        box = interval(b1[1], b1[2]) × interval(b2[1], b2[2])

        prob = calc_measure(BivDist, box)
        if prob < 0; prob = 0;end

        probInt = calc_bounds(carts, masses, box)

        bel1 = sklar_survival(a1, a2, C, box)
        if bel1 < 0; bel1 = 0; end
        bel2 = minitive_sklar(a1, a2, C, box)

        t1 =  probInt.lo <= prob
        t2 = bel1 <= prob + tolerance
        t3 = bel2 <= prob

        if any( [!t1, !t2, !t3])
            println("$box     |     $prob     |      $(probInt)      |      $(bel2)      |      $(bel1)")
            failures += 1
        end
    end
    println()
    println("Failures: $failures/$ntests")
    println()
    println()
    println()
end


#=

b1 = interval(-∞, 0.6) × interval(-∞, 0.5)
b2 = interval(-∞, 0.2) × interval(-∞, 0.3)
b3 = interval(-∞, 0.5) × interval(-∞, 0.1)
b4 = interval(0.4, 0.6) × interval(0.4, 0.6)
b5 = interval(0.4, 0.6) × interval(0.4, 0.6)
b6 = interval(0.7, 0.8) × interval(0.6, 0.9)
b7 = interval(0.5, ∞) × interval(0.5, ∞)
b8 = interval(0.1, ∞) × interval(0.1, ∞)
b9 = interval(0.7, ∞) × interval(0.8, ∞)

testSets = [b1, b2, b3, b4, b5, b6, b7, b8, b9]

#testSets = [b1, b2, b3, b4]



for box in testSets

    prob = calc_measure(BivDist, box)

    probInt = calc_bounds(carts, masses, box)

    bel1 = sklar_survival(a1, a2, C, box)
    bel2 = minitive_sklar(a1, a2, C, box)

    println("$box     |     $prob     |      $(probInt)      |      $(bel2)      |      $(bel1)")

    #@test probInt.lo <= prob
    #@test bel1 <= prob
    #@test bel2 <= prob
    #@test bel1 ≈ bel2

end

=#

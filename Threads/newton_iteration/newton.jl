#!/usr/bin/julia

if length(ARGS) != 2 && length(ARGS) != 3
  println("Two or three arguments expected")
  end

degree = nothing
nmblines = nothing
for arg in ARGS
  if arg[1] == '-'
    if arg[2] == 't'
      println("Ignoring number of threads $(arg[3:end])")
    elseif arg[2] == 'l'
      global nmblines = parse(Int64, arg[3:end])
    else
      println("Unknown argument $arg")
      exit(1)
    end
  else
    global degree = parse(Int64, arg)
  end
end

if isnothing(degree)
  println("No degree given")
  exit(1)
elseif degree <= 0
  println("Degree must be positive")
  exit(1)
end

if isnothing(nmblines)
  println("No number of lines given")
  exit(1)
elseif nmblines <= 0
  println("Number of lines must be positive")
  exit(1)
end


# In addition to the roots, we use two hypothetical attractors at 0 and
# infinity. The one at infinite has index degree, the one at zero has index
# degree+1.
const attrnmb = degree + 2;

# We compute the iteration completely (per assignment description), but we cut
# off before we write to file.
const iter_cutoff = 100

# We may precompute the exact position of all roots.
const roots = [exp(2*pi*im*dx/degree) for dx in 0:degree-1]

# The distance to roots at which we stop the iteration.
const epsilon = 1e-3

# The boundary towards infinity.
const divminre = -1e10
const divmaxre = 1e10
const divminim = -1e10
const divmaxim = 1e10

function rootiter(z::Complex{Float64})
  attr = nothing
  iter = 0
  while true
    if real(z) < divminre || real(z) > divmaxre || imag(z) < divminim || imag(z) > divmaxim
      attr = degree
      break
    end

    if abs(z) < epsilon
      attr = degree+1
      break
    end

    for (dx,root) in enumerate(roots)
      if abs(z-root) < epsilon
         attr = dx-1
         break
      end
    end
    if !isnothing(attr)
      break
    end
    
    z = z * (1. - 1. / degree) + 1. / z^(degree-1) / degree

    iter += 1
  end

  return (attr, iter < iter_cutoff ? iter : iter_cutoff-1)
end


# When trying to translate this to C, note that Julia used column-major order.
# Ignoring this yields an image that is mirrored along the diagonal.
attriters = [rootiter(zre + im*zim)
             for zre in range(-2.,2., length=nmblines),
                 zim in range(-2.,2., length=nmblines)
            ]


attrfile = open("newton_attractors_x$(degree).ppm", "w");
iterfile = open("newton_convergence_x$(degree).ppm", "w");

println(attrfile, "P3")
println(attrfile, "$nmblines $nmblines")
println(attrfile, "$attrnmb")

println(iterfile, "P3")
println(iterfile, "$nmblines $nmblines")
println(iterfile, "$iter_cutoff")

for lx in 1:nmblines
  for cx in 1:nmblines
    (attr, iter) = attriters[cx,lx]

    attr1 = mod(attr + 1*div(attrnmb,3), attrnmb)
    attr2 = mod(attr + 2*div(attrnmb,3), attrnmb)
    print(attrfile, "$attr $attr1 $attr2 ")

    print(iterfile, "$iter $iter $iter ")
  end
  println(attrfile)
  println(iterfile)
end

close(attrfile)
close(iterfile)

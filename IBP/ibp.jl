import Base: filter, print, println

type Ibp
  sets::Array{Set{Int64},1}
end

psh(s::Ibp,x::Ibp) = psh(s,x.sets)
psh(s::Ibp,x::Set{Int64}) = push!(s.sets,x)

.==(a::Ibp,b::Ibp) = compareIbp(a,b) == "="
.!=(a::Ibp,b::Ibp) = compareIbp(a,b) == "!="
.>(a::Ibp,b::Ibp) = compareIbp(a,b) == ">"
.<(a::Ibp,b::Ibp) = compareIbp(a,b) == "<"
function compareIbp(a::Ibp,b::Ibp)
  la = length(a.sets)
  lb = length(b.sets)

  if la > lb
    out = ">"
  elseif la < lb
    out = "<"
  else
    if la == 0
      out = "="
    else
      sa = (sort(a)).sets
      sb = (sort(b)).sets
      if all( {sa[i] .== sb[i] for i in 1:la} )
        out = "="
      else
        out = "!="
      end
    end
  end

  return out
end

function concat(a::Ibp,b::Ibp)
  return Ibp([a.sets,b.sets])
end

# Compares sets. Larger sets are sets with more & larger elements.
function elwiseSetCompare(a::Set,b::Set)
  la = length(a)
  lb = length(b)

  if la > lb 
    out = ">"
  elseif la < lb
    out = "<"
  else
    function compareLargest(as,bs)
      if length(as) == 0
        iout = "="
      else
        mas = maximum(as)
        mbs = maximum(bs)
        if mas > mbs 
          iout = ">"
        elseif mas < mbs
          iout = "<"
        else
          iout = compareLargest(setdiff(as,mas),setdiff(bs,mbs))
        end
      end
      return iout
    end
    out = compareLargest(a,b)
  end
  return out
end

.>(a::Set,b::Set) = elwiseSetCompare(a,b) == ">"
.<(a::Set,b::Set) = elwiseSetCompare(a,b) == "<"
.==(a::Set,b::Set) = elwiseSetCompare(a,b) == "="

# Sorts the sets within the Ibp array
sort(s::Array{Set{Int64},1},rev=false) = sort(Ibp(s),rev)
function sort(s::Ibp,rev=false)
  x = s.sets
  n = length(x)
  if n > 1
    newn = 1
    while newn > 0
      newn = 0
      for i in 2:n
        if x[i-1] .> x[i]
          temp = x[i-1]
          x[i-1] = x[i]
          x[i] = temp
          newn = i
        end
      end
      n = newn
    end # while
  end # if n > 0
  return Ibp(x)
end


function println(s::Set{Int64})
  print(s)
  println()
end
function print(s::Set{Int64})
  print("{")

  if length(s) > 0
    counter = 0
    for x in s
      counter += 1
      if counter > 1
        print(",",x)
      else 
        print(x)
      end
    end
  end

  print("}")
end

function println(a::Ibp)
  print(a)
  println()
end

function print(a::Ibp)
  s = a.sets
  la = length(s)
  print("[")
  if la > 0
    counter = 0
    for x in s
      counter += 1
      if counter > 1
        print(",")
      end
      print(x)
    end
  end
  print("]")
end

string(s::Set{Int64}) = string(Ibp([s]))
function string(a::Ibp) 
  s = Base.string(sort(a))
  s = replace(s,"Set{Int64}","")
  s = replace(s,"Ibp","")
  s = replace(s,"(","")
  s = replace(s,")","")
  s = replace(s,"[","")
  s = replace(s,"]","")
  s = replace(s,"{","")
  s = replace(s,"},","|")
  s = replace(s,"}","")
  return s
end

function tally(x)
  ux = unique(x)
  ul = length(ux)
  dict = {ux[i] => i for i in 1:ul}
  ret = zeros(Int64,ul)
  for k in x
    ret[dict[k]] += 1
  end
  return ret 
end

function tally(z::Ibp) 
  B = length(z.sets)
  u = unique(z.sets)
  tallyNums = tally(z.sets)
  pairs = {(tallyNums[i],u[i]) for i in 1:length(u)}
  pairs = Base.sort(pairs,by=x->x[1],rev=true)
  return pairs
end

function mode(z::Array{Ibp,1})
  return tally(z)[1][2]
end

function tally(z::Array{Ibp,1};prop=false)
  dict = (Array{Set{Int64},1} => Int64)[]
  for el in z
    x = sort(el).sets
    if haskey(dict,x)
      dict[x] += 1
    else 
      dict[x] = 1
    end
  end
  out = collect(dict)
  if prop
    B = length(z)
    out = { ( x[2] / B , Ibp(x[1]) ) for x in out }
  else
    out = { ( x[2] , Ibp(x[1]) ) for x in out }
  end
  sortedOut = Base.sort(out,by=x->x[1],rev=true)
  
  return sortedOut 
end

#function tally(z::Array{Ibp,1}) # ≈ 90 seconds for Ibp class; 14 seconds for 1 million Z matrices
#  println("Converting to string") # This step takes the longest
#  @time s = {string(x) for x in z}
#  println("Done converting to string")
#  u = unique(s)
#  tallyNums = tally(s)
#  pairs = {(tallyNums[i],toIbp(u[i])) for i in 1:length(u)}
#  pairs = Base.sort(pairs,by=x->x[1],rev=true)
#  return pairs
#end

function toIbp(s::String)
  as = split(s,"|")
  if as[1] == ""
    K = 0
  else
    K = length(as)
  end
  z = Array(Set{Int64},K)
  if K > 0
    for k in 1:K
      ai = split(as[k],",")
      z[k] = Set{Int64}({int(a) for a in ai})
    end
  end
  return Ibp(z)
end


function filter(f::Function,S::Ibp;col=0,rm_empty_sets=true) 
  K = length(S.sets)

  if col == 0
    sets = {filter(f,S.sets[k]) for k in 1:K}
  else
    sets = {ifelse(k==col,filter(f,S.sets[k]),S.sets[k]) for k in 1:K}
  end

  if rm_empty_sets
    sets = filter(x -> x != Set{Int64}({}), sets)
  end

  return Ibp(sets)
end

function choose(n,k)
  nCk =  exp( lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1) )
  return  convert(Int64,round(nCk))
end

function cardIbp(N,K) # counts number of possible multisets for N customers and up to K sets
  if K > 0
    S = 2^N # Number of possible sets with N customers
    card = choose(K+S-1,K) # Number of ways to assign K balls into S bins, some may be empty
  elseif K == 0
    card = 1
  end

  return card 
end

#=
  S = Ibp([Set(2,3),Set(2,3),Set(2,3)])
  S = Ibp([Set(1,2,3),Set(1,2),Set(3)])
  S = Ibp([Set(1,2,3)])
  enumNewPartition(S,2)
=#
function enumNewPartitions(S::Ibp,exclude::Int64)

  S = filter(x -> x != exclude, S)
  if length(S.sets) > 0
    ts = tally(S)
    numNewPartitions = prod({x[1] + 1 for x in ts})
    newPartitions = Array(Ibp,numNewPartitions)
    lts = length(ts)

    c = 0
    function enum(x=1,set=[])
      n = ts[x][1]
      xs = ts[x][2]
      for i in 0:n
        xset = {ifelse(i >= d, Set([collect(xs),exclude]), xs) for d in 1:n}
        if x < lts
          enum(x+1,[xset,set])
        else
          c += 1
          newPartitions[c] = Ibp([xset,set])
        end
      end
    end

    enum()
  else
    newPartitions = []
  end # if length(S.sets) > 0

  return newPartitions
end

function enumNewPartsWithNewDishes(S::Ibp,exclude::Int64;K=[0:10])
  newParts = enumNewPartitions(S,exclude)
  l = length(newParts)
  if l == 0
    l = 1
    newParts = [Ibp([])]
  end
  
  out = Array(Ibp,l*length(K))
  c = 0
  for i in 1:l
    x = newParts[i]
    for  k in K
      c += 1
      y = x
      new = [Set(exclude) for j in 1:k]
      out[c] = Ibp([y.sets, new])
    end
  end

  return out
end

function toMatrix(z::Ibp,N,colOrder=[1:length(z.sets)])
  K = length(z.sets)
  Z = zeros(Int64,N,K)

  for i in 1:N
    for k in 1:K
      Z[i,k] = ifelse(in(i,z.sets[k]),1,0)
    end
  end

  Z = Z[:,colOrder]

  return Z
end

function table(S::Ibp)
  sets = S.sets
  K = length(sets)

  if K > 0 
    dict = (Set{Int64} => Int64)[] # Empty Dictionary
    for k in 1:K
      if !haskey(dict, sets[k])
        dict[sets[k]] = 1
      else
        dict[sets[k]] += 1
      end
    end
  else
    dict = 0
  end

  return dict
end

function table(S::Ibp, i::Int64; before=true) # ONLY TO BE USED IN pold!!! Must rm new dishes!!!
  # Counts the unique sets pre-i and the the number of those sets post-i in S_{-i}
  S1 = filter(x -> x <= i, S) # S including customer i and before, assuming  i > 1
  sets = filter(x -> x != Set(i), S1.sets)  # remove {σᵢ} from S1
  K = length(sets)

  if K > 0 
    dict = ( Set{Int64} => (Int64,Int64) )[] # Empty Dictionary
    for k in 1:K
      if before
        sk0 = filter(x -> x < i, sets[k])
      else
        sk0 = filter(x -> x != i, sets[k])
      end
      sk1 = sets[k]
      if !haskey(dict, sk0)
        dict[sk0] = (1,0)
      else
        dict[sk0]= (dict[sk0][1] + 1, dict[sk0][2])
      end
      if in(i,sk1)
        dict[sk0] = (dict[sk0][1], dict[sk0][2] + 1)
      end
    end
  else
    dict = 0
  end

  return dict
end

function tableToIbp(tab::Dict{Set{Int64}, (Int64,Int64)})
  S = Ibp([])

  if tab != 0
    for t in tab
      for i in 1:t[2][1]
        psh(S, t[1])
      end
    end
  end

  return S
end

function sumMatrices(x) 
  myncol(m) = size(m,2)

  l = length(x)
  @time k=maximum(pmap(myncol,x))
  n = size(x[1],1)
  out = zeros(Int,n,k)
  for i in 1:l
    out += [x[i] zeros(Int,n,k-myncol(x[i]))]
  end
  return out
end

function lof (z) 
  lofz = z + 0
  n,k = size(z)
  if k > 1
    z = z[:,shuffle([1:k])]
    ord = {minimum(find(x->x==1,z[:,i])) for i in 1:k}
    lofz = z[:,sortperm(ord)]
  end
  return lofz
end

function meanZ (Z; p=.8)
  mz = sumMatrices(Z) / length(Z)
  mz = ifelse(mz .> p, 1, 0)
  K = size(mz,2)
  colsum = {sum(mz[:,k]) for k in 1:size(mz,2)}
  ind = colsum .> 0
  return mz[:,ind]
end

function matrixToMultiset(z) 
  N,K = size(z)
  s = Ibp([])

  for k in 1:K
    ind = find(x -> x == 1, z[:,k])
    psh( s, Set(ind) )
  end
  return s
end

#=
  include("ibp.jl")
  z = Ibp([ Set(1), Set(1), Set(1), Set(1,2), Set(2), Set(1,2), Set(1,2), Set(2) ])
  tab = table(z,2)
  tableToIbp(tab)
=#

#=
  include("ibp.jl")
  s = Set(1)
  a = Ibp([Set(1),Set(2,3,1),Set(1,4)])
  b = Ibp([Set(1),Set(2,3)])
  c = concat(a,b)
  sort(a)
  sort(b)
  sort(c) .== sort(b)
  sort(a) .== Ibp(shuffle(a.sets))
  print(a.sets[2])
  print(a)
  string(a)
=#

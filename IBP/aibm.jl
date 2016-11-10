#= useful reads for creating classes, and Sets.jl
  https://github.com/JuliaLang/julia/blob/master/base/set.jl#L25
  http://julia.readthedocs.org/en/latest/manual/constructors/
  http://quant-econ.net/jl/types_methods.html#defining-types-and-methods
  http://julia.readthedocs.org/en/latest/manual/types/
  
  Weighted Sampling:
  wsample([1,2,3],[.1,.2,.7])
  https://github.com/JuliaStats/StatsBase.jl/blob/master/src/sampling.jl
=#

println("Loading packages...")
  using Gadfly,Distributions,DataFrames
println("Done loading packages...")

include("ibp.jl") # the Ibp type: array of sets

expDecay(d) = exp(-d)

function probOld(p2::Ibp,t,f,D,sig,logout) # HERE
  pf2 = Ibp(filter(x -> x != Set(sig[t]), p2.sets)) # All sets in π₂ ( excluding {σₜ} )
  p1 = filter(x -> x != sig[t], p2) # All sets in π₁ ( excluding {σₜ} )
  cu = tally(p1)  # Unique sets in p1
  dv = tally(pf2) # Unique sets in p2 ( excluding {σₜ} )
  lu = length(cu)
  vorig = {vi[2] => filter(x -> x != sig[t],vi[2]) for vi in dv} # origin u of each unique set v
  c = {cu[i][2] => cu[i][1] for i in 1:lu} # Counts of unique sets in u indexed by unique sets in u
  d = (Set{Int64} => Int64)[] # Counts of taking new customer in next partition indexed by sets in u
  uind = {p1.sets[i] => i for i in 1:length(p1.sets)}

  for el in pf2.sets
    v = el
    u = vorig[v]
    if !haskey(d,u)
      d[u] = 0
    end
    if in(sig[t],v)
      d[u] += 1
    end
  end

  p_old = ifelse(logout,0,1)
  p = probOldDish(t,p1,f,D,sig)

  for el in cu
    u = el[2]
    if logout
      p_old += logpdf( Binomial( c[u] , p[ uind[u] ] ), d[u] )
    else
      p_old *=    pdf( Binomial( c[u] , p[ uind[u] ] ), d[u] )
    end
  end
  
  return p_old
end

function probOldDish(i,S,f,D,sig)
  lam(i,j) = f(D[i,j])
  sets = S.sets
  K = length(sets)
  m = sum( [length(s) for s in sets] )
  hnum = Array(Float64,K)
  for k in 1:K
    cust = collect(sets[k])
    hnum[k] = sum( [lam(sig[i],c) for c in cust] )
  end
  
  #h = hnum / sum(hnum)
  #p = h * m / i
  logh = log(hnum) - log( sum(hnum) )
  logp = logh + log(m) - log(i)
  p = exp(logp)
  p = min(p,1)
  
  return p
end

function ribp(;B=1, N=3, a=1, decay=expDecay)
  out = Array(Ibp, B)
  @time for b in 1:B
    out[b] = ribp(N, a, decay=decay)
  end
  
  return out
end

function ribp(N::Int64, a::Float64; decay=expDecay)
  S = Ibp([])

  for i in 1:N
    # Sample Old Dishes
    for k in 1:length(S.sets)
      p = length(S.sets[k]) / i
      if p > rand()
        push!(S.sets[k], i)
      end
    end
    x = rand(Poisson(a/i))

    # Sample New Dishes
    for xx in 1:x
      psh(S, Set(i) )
    end
  end

  return S
end

function dibp(S::Ibp, N::Int64, a; decay=expDecay, logout=false)

  function pnew (S)
    K = length(S.sets)
    x = zeros(Int64,N)

    for k in 1:K
      owner = minimum(S.sets[k])
      x[owner] += 1
    end

    p = sum({ logpdf(Poisson(a/i),x[i]) for i in 1:N })
    return p
  end

  function pold (S,i) # HERE
    function pk (set) 
      m = length(filter(x -> x < i, set))
      return m / i
    end

    U = table(S,i)
    if U != 0
      J = length(U)
      po = { logpdf( Binomial(u[2][1], pk(u[1]) ), u[2][2] ) for u in U }
    else
      po = 0
    end

    return sum(po)
  end

  po = 0
  for i in 1:N
    if i > 1
      po += pold(S,i)
    end # if i > 1
  end # for i in 1:N

  pn = pnew(S)
  pr = po + pn
  prob = ifelse(logout, pr, exp(pr) )
  return prob
end

#= 
  include("aibm.jl")
  B = 1000000; a=1.0; N = 3
  out = ribp(B=B, N=N, a=a);
  tout = tally(out,prop=true);
  {(tout[i][1],dibp(tout[i][2],N, a),tout[i][2]) for i in 1:20}
  dpois(a,x) = pdf(Poisson(a),x)
=#

function daibm(S::Ibp,a,D,decay=expDecay,sig=[1:size(D,1)],tau=1,logout=false)
  sets = S.sets
  isets = [Set{Int64}({}) for s in sets]
  N = size(D,1)
  K = length(sets)

  if K > 0
    f(d) = decay(d*tau)
    invsig = sortperm(sig) 
    x = zeros(Int64,N)
    owner = zeros(Int64,K)
    for k in 1:K
      isets[k] = Set(invsig[collect(sets[k])])
      owner[k] = sig[ minimum(isets[k]) ]
      x[owner[k]] += 1
    end
    y = cumsum(x[sig])

    p_old = ifelse(logout,0,1)
    p_new = ifelse(logout,0,1)
    for i in 1:N
      if i > 1
        if y[i-1] > 0
          Scur = filter(x -> invsig[x] <= i, S)
          pt = probOld(Scur,i,f,D,sig,logout)
          if logout
            p_old += pt
          else
            p_old *= pt
          end
        end
      end
      if logout
        p_new += logpdf(Poisson(a/i),x[sig[i]])
      else
        p_new *= pdf(Poisson(a/i),x[sig[i]])
      end
    end
    if logout
      p = p_old + p_new
    else
      p = p_old * p_new
    end

  else # i.e. if K == 0
    if logout
      p = sum({logpdf(Poisson(a/i),0) for i in 1:N})
    else
      p = prod({pdf(Poisson(a/i),0) for i in 1:N})
    end
  end

  return p
end

function raibm(a,D,decay=expDecay,sig=[1:size(D,1)],tau=1)
  f(d) = decay(d*tau)

  N = size(D,1) # number of customers
  x = {rand(Poisson(a/i)) for i in 1:N} # number of new dishes taken
  K = sum(x)
  y = cumsum(x)
  #S = Ibp( [Set{Int64}({}) for k in 1:K] )
  S = Ibp( [] )

  if K > 0
    for i in 1:N
      # sample old dishes
      if i > 1 && y[i-1] > 0
        p = probOldDish(i,S,f,D,sig)
        for k in 1:y[i-1]
          if p[k] > rand()
            push!(S.sets[k],sig[i])
          end
        end
      end

      # sample new dishes
      if x[i] > 0
        for k in 1:x[i]
          psh(S,Set(sig[i]))
        end
      end # if x[i] > 0
    end # for i in 1:N
  end

  return S
end

# 16 July, 2015
function gibbs(B,a,D,decay=expDecay,sig=[1:size(D,1)],tau=1)
  isig = invperm(sig)
  N = length(sig)
  z = Array(Ibp,B)
  z[1] = Ibp([])
  #z[1] = Ibp([Set(1), Set(2), Set(3,1)])

  for b in 2:B
    z[b] = z[b-1]
    for i in 1:N
      #parts = enumNewPartitions(z[b],sig[i])
      parts = enumNewPartsWithNewDishes(z[b],sig[i],K=(0:5))
      lparts = length(parts)

      probs = zeros(lparts)
      for lp in 1:lparts
        probs[lp] = exp( daibm(parts[lp],a,D,decay,sig,tau,true) )
      end

      ind = wsample(1:lparts, probs)
      z[b] = parts[ind]
    end # i in 1:N

    if b%(B/1000)==0 print("\r",round(100*b/B,1),"%") end
  end # b in 2:B

  return z
end

function gibbs2(B,a,D;decay=expDecay,sig=[1:size(D,1)],tau=1)
  N = length(sig)
  z = Array(Ibp,B)
  z[1] = Ibp([])
  
  for b in 2:B
    zz = Ibp( collect(z[b-1].sets) ) 
   
    for i in 1:N
      K = length(zz.sets)

      # Sample old dishes
      for k in 1:K
        mk = length(filter(x -> x != sig[i], zz.sets[k]))

        z0 = filter(x -> x != sig[i], zz, col=k, rm_empty_sets=false)
        z1 = filter(x -> x > 0, z0, rm_empty_sets=false)
        z1.sets[k] = Set([ collect(z1.sets[k]), sig[i] ])

        fz0 = filter(x -> x > 0, z0)
        fz1 = filter(x -> x > 0, z1)

        c = sum(z0.sets[k] .== z0.sets)
        d = sum(z1.sets[k] .== z1.sets)

        lp0 = daibm(fz0,a,D,decay,sig,tau,true) + log(c)
        lp1 = daibm(fz1,a,D,decay,sig,tau,true) + log(d)

        mx = max(lp0, lp1)
        p0 = exp(lp0 - mx)
        p1 = ifelse(mk == 0, 0, exp(lp1 - mx)) #new

        ind = wsample([0,1], [p0,p1])
        zz = ifelse(ind == 1, z1, z0)
      end
      zz = filter(x -> x > 0, zz)

      # Sample new dishes
      #x = rand(Poisson(a/N))
      #for k in 1:x
      #  psh( zz, Set(sig[i]) )
      #end

      nd = 3
      x = [0:nd]
      S = { Ibp([ zz.sets, [] ]) for xx in x }
      pNew = zeros(length(x))
      ind = 0
      for xx in x
        ind += 1
        for k in 1:xx
          psh( S[ind], Set(sig[i]) )
        end
        pNew[ind] = logpdf(Poisson(a/N),xx)#daibm(S[ind],a,D,decay,sig,tau,true)
      end
      mx = maximum(pNew)
      pNew = pNew - mx
      #println(S[1])
      #println(S[2])
      #println(S[3])
      #println(S[4])
      #println( exp(pNew) / sum(exp(pNew)) )
      #readline()
      p = wsample( x , exp(pNew) )
      zz = Ibp([ S[p+1].sets, [] ])

    end # i in 1:N
    z[b] = zz

    if b%(B/1000)==0 print("\r",round(100*b/B,1),"%") end
  end # b in 2:B

  return z
end

#= 
  include("aibm.jl")
  #D = [[0 1 1]; [1 0 1]; [1 1 0]]; B = 10000; a = 1; sig = [3,2,1]
  D = [[0 1 2]; [1 0 3]; [2 3 0]]; B = 10000; a = 1; sig = [3,2,1]; tau = 3
  @time out = gibbs2(B,a,D,sig=sig,tau=tau);
  mean({length(x.sets) for x in out})

  @time tout = tally(out,prop=true);
  {(tout[i][1], round(daibm(tout[i][2],a,D,expDecay,sig,tau),4), tout[i][2]) for i in 1:20}
  sum({daibm(t[2],a,D,expDecay,sig,tau) for t in tout})
=#

function gibbs3(B,N,a,decay=expDecay) #HERE
  z = Array(Ibp,B)
  z[1] = Ibp([])
  
  for b in 2:B
    z[b] = zz = Ibp(collect(z[b-1].sets))
    x = rand(Poisson(a/N), N)

    for i in 1:N
      # Sample old dishes
      K = length(zz.sets)
      for k in 1:K
        mk = length(filter(x -> x != i, zz.sets[k]))
        #pk = mk / N

        z0 = filter(x -> x!= i, zz, col=k, rm_empty_sets=false)
        z1 = Ibp([ collect(zz.sets) ])
        z1.sets[k] = Set([collect(z1.sets[k]),i])

        f0 = filter(x -> x > 0, z0)
        f1 = filter(x -> x > 0, z1)

        c = sum(z0.sets[k] .== z0.sets)
        d = sum(z1.sets[k] .== z1.sets)
        lp0 = dibp(f0,N,a,logout=true) + log(c)
        lp1 = dibp(f1,N,a,logout=true) + log(d)

        mx = max(lp0, lp1)
        p0 = exp(lp0 - mx)
        p1 = exp(lp1 - mx)

        pkk = p1 / (p1+p0)
        pkk = ifelse(mk == 0, 0, pkk)
        
        #if round(pkk,5) != round(pk,5)
        #  print("z0: "); print("p0: ", exp(lp0)); println(z0)
        #  print("z1: "); print("p1: ", exp(lp1)); println(z1)
        #  println("pk: ",pk,"; ","pkk: ",pkk)
        #  readline()
        #end

        if pkk > rand()
          #zz.sets[k] = Set([collect(zz.sets[k]),i])
          zz = z1
        else
          #zz = filter(x -> x != i, zz, col=k, rm_empty_sets=false)
          zz = z0
        end
      end

      # Sample new dishes
      zz = filter(x -> x != 0, zz)
      for k in 1:x[i]
        psh(zz, Set(i))
      end
    end # i in 1:N

    z[b] = Ibp(collect(zz.sets))

    if b%(B/1000)==0 print("\r",round(100*b/B,1),"%") end
  end # b in 2:B

  return z
end

#=
  include("aibm.jl")
  B = 100000; N = 3; a = 1
  @time out = gibbs3(B,N,a);
  @time tout = tally(out,prop=true);
  {(tout[i][1], dibp(tout[i][2],N,a), tout[i][2]) for i in 1:10}
  sum({dibp(t[2],N,a) for t in tout})
=#

function mh (a,D,decay=expDecay,sig=[1:size(D,1)],tau=1;B=1000000)
  
  isig = invperm(sig)
  N = length(sig)
  z = Array(Ibp,B)
  z[1] = Ibp([])
  acc = 0
  counter = 0

  for b in 2:B
    S = z[b-1]
    for i in 1:N
      S_cand = filter(x -> x != sig[i], S)

      # Propose existing dishes
      K = length(S_cand.sets)
      if K > 0
        for k in K
          if .5 > rand()
            push!(S_cand.sets[k],sig[i])
          end
        end
      end

      # Prospose new dishes
      #x = rand(Poisson(a/i))
      x = sample([0:(floor(a)+1)])
      if x > 0
        for k in 1:x
          psh(S_cand,Set(sig[i]))
        end
      end
      
      mhratio = daibm(S_cand,a,D,decay,sig,tau,true) - daibm(S,a,D,decay,sig,tau,true)
      if mhratio > log(rand())
        S = S_cand
        acc += 1
      end
    end # end of i in 1:N

    z[b] = S

    if b%(B/1000)==0 print("\r",round(100*b/B,1),"%") end
  end # end of b in 1:B 

  println("Acceptance Rate: ", round(100 * acc / (N*B),3), "%")
  return z
end

#=
  include("aibm.jl")
  a = 1
  D = [[0 1 3]; [1 0 2]; [3 2 0]]

  # USArrest
  D = [[0.00 0.12 0.66 3.74 3.78];
       [0.12 0.00 0.59 3.65 3.69];
       [0.66 0.59 0.00 3.32 3.46];
       [3.74 3.65 3.32 0.00 1.01];
       [3.78 3.69 3.46 1.01 0.00]]

  #sig = [3,1,2]
  sig = [1:size(D,1)]
  tau = 0
  
  B = 100000
  z = Array(Ibp,B)
  @time for b in 1:B # ≈ 30 seconds for a million, 70 seconds with sort.
    z[b] = raibm(a,D,expDecay,sig,tau)
    if b%(B/10)==0 print("\r",round(100*b/B),"%") end
  end

  @time tz = tally(z) # a million draws: 30 seconds
  {(tz[i][1]/B, daibm(tz[i][2],a,D), tz[i][2]) for i in 1:20}

  sigDiff(vals) = abs((vals[1] - vals[2]) / sqrt(vals[1]*(1-vals[1])/B)) > 1.96

  bb = length(tz)
  sigdiff = zeros(bb)
  @time for ind in 1:bb
    vals = ( tz[ind][1] / B, exp(daibm(tz[ind][2],a,D,expDecay,sig,tau,true)) )
    sigdiff[ind] = sigDiff(vals)
    if ind%(div(bb,10))==0 print("\r",round(100*ind/bb),"%") end
  end
  mean(sigdiff)

  emp = zeros(bb)
  @time for ind in 1:bb
    exp(daibm(tz[ind][2],a,D,expDecay,sig,tau,true))
    if ind % 10 == 0 print("\r",round(100*ind/bb),"%") end
  end
  
  bad = find(x -> x == 1,sigdiff)
  for ind in bad
    print(  "Partition:   "); println(tz[ind][2])
    println("Empirical:   ",  tz[ind][1]/B)
    p = exp(daibm(tz[ind][2],a,D,expDecay,sig,tau,true))
    println("Theoretical: ",  p)
    println("Ratio:       ",  tz[ind][1]/B / p)
    println()
  end

  for i in bad
    print(  "Partition",i,":   "); println(tz[i][2])
    println(tz[i][1]/B, "\t ", exp( daibm(tz[i][2],a,D,expDecay,sig,tau,true) ))
    println()
  end

  include("aibm.jl")
  a = 10
  D = [[0 1 1]; [1 0 1]; [1 1 0]]
  sig = [1:size(D,1)]
  tau = 1
  B = 10000
  @time out = mh(a,D,expDecay,sig,tau,B=B)
  @time tal = tally(out)
  @time { (tal[i][1]/B, daibm(tal[i][2],a,D,expDecay,sig,tau), tal[i][2]) for i in 1:length(tal)}
  sum({daibm(tal[i][2],a,D,expDecay,sig,tau) for i in 1:length(tal)})

  z = Array(Ibp,B)
  @time for b in 1:B # ≈ 30 seconds for a million, 70 seconds with sort.
    z[b] = raibm(a,D,expDecay,sig,tau)
    if b%(B/10)==0 print("\r",round(100*b/B),"%") end
  end
  @time tal2 = tally(z)
  sum({daibm(tal2[i][2],a,D,expDecay,sig,tau) for i in 1:length(tal2)})

  umh = {t[2] for t in tal}
  uz = {t[2] for t in tal2}

  mean({length(u.sets) for u in umh})
  mean({length(u.sets) for u in uz})

=#

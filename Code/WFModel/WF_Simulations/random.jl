using Distributed
@everywhere using Distributions
@everywhere using Random
using NPZ
using Plots
using SpecialFunctions


@everywhere function p_simula(N,fit,m)   #one simulation for one landscape
    sel = zeros(length(fit))
    pop = zeros(length(fit))
    pop[1] = N  
    while true
        
        #MUTATION
        pop .+= pop[1] .* m
        pop[1] += N - sum(pop) 
        
        #SELECTION
        sel .= fit .* pop       
        sel ./= sum(sel)
                        
        #GEN DRIFT
        dis = Multinomial(N, sel)
        rand!(dis,pop)
        
        if pop[1] <= 0.1*N && maximum(pop) > 0.8*N
            break
        end 
    end
    findmax(pop)[2]
end

function mutation_m(mus,mud,ns,nd)   #creating mutation matrix
    m = zeros(ns + nd + 1)
    for i in 2:ns+1
        m[i] = mus
    end
    for i in ns+2:ns+nd +1
        m[i] = mud
    end
    #m[1] = 1 - (ns*mus + nd*mud)
    m
end

@everywhere function simula_tot(RUN,N,fit,m)  #average for one landscape
    ris = zeros(RUN)
    for i in 1:RUN
        ris[i] = p_simula(N,fit,m)
    end  
    gt = 1:length(fit)
    freq = [sum(ris .== b) for b in gt]
    freq /= sum(freq)
    freq
end

if length(ARGS) > 0
    seed = parse(Int64,ARGS[1])
else
    seed = rand(1:10000000)
end
Random.seed!(seed)


alpha = 4*10^(-3)
ns = 12                                      # benficial singles 
nd = 20                                      # benficial doubles 

mus = 10^(-10)*3*0.76/5.80              # SN mutation rate (locus mutation rate)
mud = 10^(-10)*alpha*3*0.99*0.52/9.63      # DN mutation rate (locus pair mutation rate)


U = ns*mus+nd*mud

#meanExp = 0.1
#f = rand(Exponential(meanExp),1+ns+nd) .+ 1.   # fitness values

#meanExpS = 0.1
#meanExpD = 0.1
#fS = rand(Exponential(meanExpS), 1+ns) .+ 1.
#fD = rand(Exponential(meanExpD), nd) .+ 1.
#f = vcat(fS, fD)

#meanWeibull = 0.01
#shape = 1.0
#scale = meanWeibull/gamma(1+1/shape)
#p = Weibull(shape,scale)
#f = rand(p, 1+ns+nd) .+ 1.


meanParetoS = 0.01
meanParetoD = 0.01
shape = 1.5
scaleS = (meanParetoS*(shape-1))/shape
scaleD = (meanParetoD*(shape-1))/shape
fS = rand(Pareto(shape, scaleS), 1+ns) .+ 1.
fD = rand(Pareto(shape, scaleD), nd) .+ 1.
f = vcat(fS, fD)
#p = Pareto(shape,scale)
#f = rand(p, 1+ns+nd) .+ 1.


#p = PGeneralizedGaussian(0, 1960000, 0.3);
#p = PGeneralizedGaussian(0, 100, 1.0)
#x = rand(p, 1+ns+nd)
#f = broadcast(abs,x).+ 1.




f[1] = 1.                                    # fitness WT
f0 = copy(f[1:1+ns])

#f = ones(1+ns+nd)
#f[ns+2:nd+ns+1] .+= rand(Exponential(0.1),nd)
#f[2:ns+1] .= f[ns+2:2*ns+1]

findmax(f)

m = mutation_m(mus,mud,ns,nd)                # creating mutation matrix
m0 = ones(ns+1)*mus
m0[1] = 0.

#---------------- SIMULATION ------------------
#compilation
n = trunc(Int,2/U)                                 # 
p = nworkers()                                     # 
_ = simula_tot(2,n,f,m)                            # 
pmap(x -> simula_tot(x,n,f,m), [2 for i in 1:p])   # 

#run
RUN = 1000_000                                     # Averages per landscape
run = trunc(Int, RUN/p)                            # 

MS = [10^i for i in range(-3.5,2,60)]                # Log Mutation supply
Y = zeros(length(MS))                              # Y-axis (F_DN)
data = zeros(length(MS),length(f))                 # Y-axis (F_DN) for each landscape
for i in 1:length(MS)                              
    ms = MS[i]
    N = trunc(Int,ms/U)
    ris = pmap(x -> simula_tot(x,N,f,m), [run for i in 1:p])

    ris = mapreduce(permutedims, vcat, ris)
    ris = mean(ris, dims=1)
    data[i,:] = ris
    Y[i] = sum(ris[2+ns:1+ns+nd])
end

data0 = zeros(length(MS),length(f0))
for i in 1:length(MS)
    ms = MS[i]
    N = trunc(Int,ms/U)
    ris = pmap(x -> simula_tot(x,N,f0,m0), [run for i in 1:p])
    ris = mapreduce(permutedims, vcat, ris)
    ris = mean(ris, dims=1)
    data0[i,:] = ris
 
end


#scatter(MS,Y,xscale=:log10,lw = 3,label=false, color="blue", dpi = 300)
#plot!(MS,Y,xscale=:log10,lw = 3,label=false, color="blue")

#---------------- SAVING ------------------
path = "dataRandom/$seed"
npzwrite(path*"_ris.npy", data)
npzwrite(path*"_ris0.npy", data0)
npzwrite(path*"_MS.npy", MS)
npzwrite(path*"_fit.npy", f)
#savefig(path*"_plot.png")


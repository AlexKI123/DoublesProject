using Distributions
using Random
using NPZ
#using Plots

#using ProgressBars

function p_simula(N,fit,m)
    sel = zeros(length(fit))
    pop = zeros(length(fit))
    pop[1] = N  
    cont = 0
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
        if cont > 10^5
            break
        end
        cont += 1
    end
    findmax(pop)[2]
end

function mutation_m(mus,mud,ns,nd)
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

function simula_tot(RUN,N,fit,m)
    ris = zeros(RUN)
    Threads.@threads for i in 1:RUN
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


#fs = npzread("singles_TEM1.npy")
#fd = npzread("doubles_TEM1.npy")

fs = npzread("singles_P53.npy")
fd = npzread("doubles_P53.npy")


f = vcat(1,fs,fd)
f0 = vcat(1,fs)


alpha = 4*10^(-3)
#alpha = 10^(-2)
ns = length(fs)
nd = length(fd)
#mus = 10^(-10)/3
#mud = 10^(-10)*α/9

mus = 10^(-10)*3*0.76/5.80
mud = 10^(-10)*alpha*3*0.99*0.52/9.63

#mus = 8.12*10^(-6)*3*3*0.76/5.80
#mud = 8.12*10^(-6)*3*alpha*3*0.99*0.52/9.63


U = ns*mus + nd*mud

m = mutation_m(mus,mud,ns,nd)
m0 = ones(ns+1)*mus
m0[1] = 0.

#---------------- SIMULATION ------------------

#run
RUN = 5000

MS = [10^i for i in range(-3,2.2,40)]

data = zeros(length(MS),length(f))
data0 = zeros(length(MS),length(f0))

y = zeros(length(MS))
for i in 1:length(MS)
    ms = MS[i]
    N = trunc(Int,ms/U)
    ris = simula_tot(RUN,N,f,m)
    data[i,:] = ris

    y[i] = sum(ris[ns+2:length(ris)])    
end

for i in 1:length(MS)
    ms = MS[i]
    N = trunc(Int,ms/U)
    ris = simula_tot(RUN,N,f0,m0)
    data0[i,:] = ris
 
end

#print("Plotting")
#scatter(MS,y,xscale=:log10)
#plot(MS,y,xscale=:log10)
#savefig("test.png")


#---------------- SAVING ------------------
path = "dataP53/$seed"
npzwrite(path*"_ris.npy", data)
npzwrite(path*"_ris0.npy", data0)
npzwrite(path*"_MS.npy", MS)
npzwrite(path*"_fit.npy", f)
#savefig(path*"_plot.png")


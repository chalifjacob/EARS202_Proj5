## Project 4
## EARS 202
## Jacob Chalif
using Pkg
Pkg.activate(".")
Pkg.add(["FFTW","Plots","StatGeochem"])
using FFTW, Plots, StatGeochem


## Part 1
data = importdataset("BarHarborWL2013.csv", ',', importas=:Tuple)
# Signal: water Level
wl = data.Water_Level
# Number of sample points
N = length(data.Hour_Count)
# Sample period
Ts = 1/24
# Start and end time
t0 = 0
tmax = t0 + N * Ts
# Time coordinate
t = t0:Ts:tmax
# Fast Fourier Transform
ft = fft(wl)
F = fftshift(ft)
freqs = fftshift(fftfreq(N, 1.0/Ts))
# Plot results
time_domain = plot(t[1:end-1], wl, 
    title="Signal", 
    framestyle=:box,
    label = "",
    xlabel = "Days",
    ylabel = "Water level")
time_domain_cropped = plot(t[1:end-1], wl, 
    title="Signal, zoomed", 
    framestyle=:box,
    label = "",
    xlims = (5,15),
    xlabel = "Days")
freq_domain= plot(freqs, abs.(F), 
    title="Spectrum", 
    yscale=:log10, 
    framestyle=:box,
    label = "",
    xminorticks=5,
    xlabel = "Days",
    ylabel = "Spectral power")
    
# Plots in one figure
l = @layout [a b; c] # defines two-panel layout
pFin=plot(time_domain, time_domain_cropped, freq_domain, layout = l, dpi = 500, size=(600,450))
savefig(pFin, "Figures/1.png") # saves figure



## Part 2
data = importdataset("USH00273850.FLs.52j.tavg", '\t', importas=:Tuple)





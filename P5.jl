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
freq_domain= plot(1 ./ freqs, abs.(F), 
    title="Spectrum", 
    yscale=:log10, 
    framestyle=:box,
    label = "",
    xlabel = "Days",
    ylabel = "Spectral power",
    xlims = (0,5),
    xticks = 0:1:5)
    
# Plots in one figure
l = @layout [a b; c] # defines two-panel layout
pFin=plot(time_domain, time_domain_cropped, freq_domain, layout = l, dpi = 500, size=(600,450))
savefig(pFin, "Figures/1.png") # saves figure



## Part 2
# Using ERA5 6-hourly temperature data for the grid cell containing Hanover_ERA5
# Temperature is in °C
# Data spans 1940 to 2022
# I extracted this grid cell from a larger dataset using scripts I have in Matlab
data = importdataset("Hanover_ERA5.csv", ',', importas=:Tuple)
T = data.T # temperature
# Number of sample points
N = length(T)-1
# Sample period
Ts = 6/24
# Start and end time
t0 = 0
tmax = t0 + N * Ts

# time coordinate
t = t0:Ts:tmax

# Fast Fourier Transform
F = fftshift(fft(T))
freqs =  fftshift(fftfreq(length(t), 1/Ts))

# Plot results
time_domain = plot(t, T, 
    title="Signal", 
    framestyle=:box,
    label = "",
    xlabel = "Days after 1/1/40",
    ylabel = "Temperature (°C)")
time_domain_cropped = plot(t, T, 
    title="2 years", 
    framestyle=:box,
    label = "",
    xlims = (0,365*2),
    xlabel = "Days after 1/1/40")
time_domain_cropped2 = plot(t, T, 
    title="1 month", 
    framestyle=:box,
    label = "",
    xlims = (0,30),
    xlabel = "Days after 1/1/40")
freq_domain= plot(1 ./ freqs, abs.(F), 
    title="Spectrum", 
    yscale=:log10, 
    framestyle=:box,
    label = "",
    xminorticks=5,
    xlabel = "Days",
    ylabel = "Spectral power",
    xlims = (0,500))
freq_domain_cropped= plot(1 ./ freqs, abs.(F), 
    title="Spectrum, 3 days", 
    yscale=:log10, 
    framestyle=:box,
    label = "",
    xticks = -3:1:3,
    xminorticks=2,
    xlabel = "Days",
    ylabel = "Spectral power",
    xlims = (0,3))
    
# Plots in one figure
l = @layout [a b c; d e] # defines two-panel layout
pFin=plot(time_domain, time_domain_cropped, time_domain_cropped2, freq_domain, freq_domain_cropped,layout = l, dpi = 500, size=(600,450))
savefig(pFin, "Figures/2.png") # saves figure
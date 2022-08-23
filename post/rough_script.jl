using CairoMakie, DelimitedFiles

simulation_name = "test_01/"
folder_name = "seas-seme-devel/"
out_path = joinpath(dirname(pwd()), folder_name, "data/", simulation_name)
mesh_path = joinpath(out_path, "mesh/")
fig_path = joinpath(dirname(pwd()), folder_name, "plots/", simulation_name)

mkpath(fig_path)
const dpi = 300
my_theme = Theme(
    fontsize=48,
    labelsize=36,
    linewidth=6
)
set_theme!(my_theme)

time_series_1d = readdlm(joinpath(out_path, "time_series_1d.out"), skipstart=1);

yr2sec = 365*24*60*60
t = time_series_1d[:,1]/yr2sec
vmax = time_series_1d[:,2]

function sliprate_plot(vmax, t)
    #%%
    figsize = (8.2, 3.35) # in inches
	f = Figure(resolution=figsize.*dpi)
	ax = Axis(f[1,1],
			title="Max. sliprate on fault",
			xlabel="Time (years)",
			ylabel="Sliprate (m/s)",
			yscale=log10)
	
	lines!(ax, t, vmax)
	save(joinpath(fig_path, "max_sliprate.png"), f)
	f
    #%%
end
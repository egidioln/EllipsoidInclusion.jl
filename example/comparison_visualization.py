import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

inside = True
good_res = 1 if inside == True else 0

type = "in" if inside else "out"
data = pd.read_csv(f"data_{type}.csv")
data["n"] = data["n"].astype('int')
res = ["our_res","sdpa_res","mosek_res"]
time = ["our_time","sdpa_time","mosek_time"]
mem = ["our_mem","sdpa_mem","mosek_mem"]

plot_name = {"our_time":"Algorithm~1","sdpa_time":"LMIs+SDPA","mosek_time":"LMIs+Mosek"}
handlers=data.boxplot(column=time, by="n", layout=(1,3),figsize=[5,2.4],return_type='dict')
# for ax, t in zip(axs,time):

# 	ax.set_title("")

ax = plt.subplot(1,3,1)
ax.set_ylabel("Execution time (s)")
for (i,t) in enumerate(time):
	ax = plt.subplot(1,3,1+i)
	ax.set_title(plot_name[t])	
	ax.set_yscale('log')
	ax.set_xlabel("$n$")
	ax.grid(linestyle="-", linewidth=0.1)
for h in handlers:
	for f in h["fliers"]:
		f.set_markersize(3)
	for m in h["medians"]:
		m.set_linewidth(2)
		m.set_color("#B04030")
	for b in [*h["boxes"], *h["whiskers"]]:
		b.set_color("#202090")


plt.suptitle("$\\mathcal{E} \\subset {\\rm int}(\\mathcal{E}_0)$" if inside else "$\\mathcal{E} \\not\\subseteq \\mathcal{E}_0 $")
plt.tight_layout()
plt.savefig(f'boxplot_{type}.eps', format='eps')

n = data["n"].unique()
time_data = data.groupby("n")[time].mean().to_numpy()
#   aux
# [x;ln(c)] =[ln(n); 1]# *ln y 
#aux = np.linalg.pinv(np.mat([n,np.ones(len(n))])).transpose()@np.log(time_data)

print(data[(data!=good_res)].groupby("n").count()[res]) # % of failure 
print(data.groupby("n").mean()[mem]/1024) # KB allocated

print(time_data[:,1:3].sum(axis=1)/time_data[:,0]/2-1) # x times faster (in avg)

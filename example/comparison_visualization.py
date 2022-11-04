import pandas as pd
import matplotlib.pyplot as plt



data = pd.read_csv("data_out.csv")

axs=data.boxplot(column=["our_time","sdpa_time","mosek_time"], by="n", layout=(1,3),figsize=[6,4])
for ax in axs:
  ax.set_yscale('log')
plt.tight_layout()
plt.show()
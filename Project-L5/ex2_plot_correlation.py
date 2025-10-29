import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("data/gc_vs_time.csv")

print(data)
corr = data["gc_percent"].corr(data["avg_time_seconds"])
print(f"\nCorrelation coefficient (GC% vs time): {corr:.3f}")

plt.scatter(data["gc_percent"], data["avg_time_seconds"], color="blue")
plt.xlabel("GC% of Sequence")
plt.ylabel("Average Reconstruction Time (s)")
plt.title(f"GC% vs Reconstruction Time (r = {corr:.3f})")
plt.grid(True)
plt.show()

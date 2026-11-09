
from pathlib import Path
import matplotlib.pyplot as plt

# 1. Locate the E. coli genome file
genome_file = Path(
    "data/GCF_000005845.2_ASM584v2_genomic.fna"
)

# 2. Read the FASTA sequence
sequence = ""

with open(genome_file, "r") as file:
    for line in file:
        if not line.startswith(">"):
            sequence += line.strip().upper()

# 3. Calculate cumulative GC skew
g_count = 0
c_count = 0
cumulative_skew = 0

minimum = float("inf")
maximum = float("-inf")

min_position = 0
max_position = 0

positions = []
values = []

for position, base in enumerate(sequence, start=1):

    if base == "G":
        g_count += 1
        cumulative_skew += 1

    elif base == "C":
        c_count += 1
        cumulative_skew -= 1

    if cumulative_skew < minimum:
        minimum = cumulative_skew
        min_position = position

    if cumulative_skew > maximum:
        maximum = cumulative_skew
        max_position = position

    # Save one point every 1000 bases
    if position % 1000 == 0:
        positions.append(position)
        values.append(cumulative_skew)

# 4. Display results
print("E. coli GC Skew Analysis")
print("-----------------------")
print("Genome length:", len(sequence))
print("G count:", g_count)
print("C count:", c_count)
print("Minimum skew:", minimum)
print("Minimum position:", min_position)
print("Maximum skew:", maximum)
print("Maximum position:", max_position)

# 5. Plot the graph
plt.figure(figsize=(12, 6))

# Plot cumulative GC skew
plt.plot(positions, values, label="Cumulative GC Skew")

# Highlight the known OriC region
plt.axvspan(
    3923767, 3924025,
    color="green",
    alpha=0.3,
    label="Known OriC Region"
)

# Mark our predicted OriC position
plt.axvline(
    3925597,
    color="red",
    linestyle="--",
    label="Predicted OriC Position"
)

plt.xlabel("Chromosome Position (bp)")
plt.ylabel("Cumulative GC Skew")
plt.title("E. coli Cumulative GC Skew Analysis")

plt.legend()
plt.grid(True)

plt.savefig("results_gc_skew.png", dpi=300)
plt.show()

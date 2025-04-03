import os
import matplotlib.pyplot as plt

# Specify the exact path of the .histo file
histo_file = os.path.join("/media/Raid/Wee/WeeYeZhi/output/Jellyfishresults/jellyfish_k_21", "k_21_mer_counts.histo")

# Check if the file exists
if not os.path.exists(histo_file):
    raise FileNotFoundError(f"Error: The file '{histo_file}' does not exist.")

def make_plot(coverage, frequency, cov_peak):
    plt.plot(coverage, frequency)
    plt.xlabel("Coverage")
    plt.ylabel("Frequency")
    plt.axvline(x=cov_peak, color='r', linestyle="--", label="Coverage Peak")
    plt.legend()
    plt.savefig("genome.png")
    plt.show()

def detect_start(frequency):
    start = 0
    for idx in range(len(frequency)):
        test_list = frequency[idx:idx + 10]
        if all(i < j for i, j in zip(test_list, test_list[1:])):
            start = idx
            break
    print(f"Start at: {start}")
    return start

def detect_end(frequency):
    end = -1
    for idx in range(len(frequency) - 1, 9, -1):
        test_list = frequency[idx - 10:idx]
        if all(i > j for i, j in zip(test_list, test_list[1:])):
            end = idx
            break
    print(f"End at: {end}")
    return end

def estimate_coverage_peak(coverage, frequency):
    higher_freq = max(frequency)
    coverage_peak = coverage[frequency.index(higher_freq)]
    print(f"Coverage Peak: {coverage_peak}")
    return coverage_peak

def estimate_genome_size(coverage, frequency, coverage_peak):
    area_under_the_curve = sum([c * f for c, f in zip(coverage, frequency)])
    genome_size = int(area_under_the_curve / coverage_peak)
    print(f"Estimated Genome Size: {genome_size} bp")
    return genome_size

def main(dataset):
    coverage = [int(entry[0]) for entry in dataset]
    frequency = [int(entry[1]) for entry in dataset]

    start = detect_start(frequency)
    end = detect_end(frequency)

    coverage = coverage[start:end]
    frequency = frequency[start:end]

    coverage_peak = estimate_coverage_peak(coverage, frequency)
    genome_size = estimate_genome_size(coverage, frequency, coverage_peak)

    make_plot(coverage, frequency, coverage_peak)

if __name__ == "__main__":
    with open(histo_file) as fh:
        data = fh.readlines()

    dataset = [entry.strip().split(" ") for entry in data]
    main(dataset)
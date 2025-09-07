#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
class SyntheticDataGenerator:
    def __init__(self, num_samples, mean_range):
        self.num_samples = num_samples
        self.mean_range = mean_range
    
    def generate_data(self):
        mean1 = np.random.uniform(*self.mean_range)
        std1 = np.random.uniform(0.4, 1.6)
        data1 = np.random.normal(mean1, std1, self.num_samples)
        mean2 = np.mean(data1)
        std2 = np.std(data1)
        print("Mean:", mean2)
        print("Standard Deviation:", std2)
        return data1
    
    def plot_distributions(self):
        data = self.generate_data()
        plt.hist(data, bins=30, alpha=0.5, color='blue')
        plt.xlabel('Value')
        plt.ylabel('Frequency')
        plt.title('Generated Data Distribution')
        plt.show()
        
    def write_to_file(self, filename):
        data = self.generate_data()
        np.savetxt(filename, data)

if __name__ == "__main__":
    generator = SyntheticDataGenerator(num_samples=80, mean_range=(3, 6))
    generator.plot_distributions()
    generator.write_to_file('synthetic_data.txt')
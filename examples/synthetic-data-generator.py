#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
class SyntheticDataGenerator:
    def __init__(self, num_samples, mean_range):
        self.num_samples = num_samples
        self.mean_range = mean_range
    
    def generate_data(self):
        mean1 = np.random.uniform(*self.mean_range)
        mean2 = np.random.uniform(*self.mean_range)
        std1 = np.random.uniform(0.5, 1.5)
        std2 = np.random.uniform(0.5, 1.5)
        data1 = np.random.normal(mean1, std1, self.num_samples)
        data2 = np.random.normal(mean2, std2, self.num_samples)
        return np.concatenate((data1, data2))
    
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
    generator = SyntheticDataGenerator(num_samples=20, mean_range=(3, 20))
    generator.plot_distributions()
    generator.write_to_file('synthetic_data.txt')

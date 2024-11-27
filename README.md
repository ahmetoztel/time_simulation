# time_simulation
# Persistent Homology and Betti Number Computation

This Python repository contains code for calculating Betti numbers and Euler characteristics based on persistent homology using simplicial complexes. The program processes random binary matrices to simulate and analyze the relationship between black pixel counts and computational performance.

## Features

- **Neighbor Detection**: Identifies neighboring black pixels in a binary matrix.
- **Simplicial Complex Construction**: Generates 0-simplices, 1-simplices, 2-simplices, and 3-simplices.
- **Boundary Matrix Generation**: Constructs sparse boundary matrices for simplices.
- **Smith Normal Form Approximation**: Approximates Smith Normal Form for matrix analysis.
- **Betti Numbers Calculation**: Computes Betti numbers and Euler characteristics.
- **Simulation and Analysis**: Analyzes execution time versus black pixel counts with a log-log plot.

## Requirements

To run the code, install the following Python libraries:

- `numpy`
- `matplotlib`
- `scipy`
- `sympy`

You can install these dependencies using pip:

```bash
pip install numpy matplotlib scipy sympy
Usage
Clone the repository:

bash
Kodu kopyala
git clone https://github.com/yourusername/your-repo-name.git
cd your-repo-name
Run the Python script:

bash
Kodu kopyala
python script_name.py
Specify the number of iterations when prompted:

bash
Kodu kopyala
Enter the number of iterations: 10
Simulation Details
Matrix Size: 100x100
Black Pixel Range: Between 1000 and 8000 black pixels per matrix.
Outputs:
Execution time vs. black pixel count graph.
Total simulation runtime displayed in seconds.
Example Output
The program generates a log-log plot displaying the execution time as a function of the number of black pixels in the binary matrix.

Sample Plot
(A sample image of the log-log plot)

Code Overview
Core Functions
find_black_pixel_neighbors(binary_matrix) Detects neighbors of black pixels.

generate_simplices(binary_matrix) Constructs simplicial complexes (0, 1, 2, and 3 simplices).

generate_boundary_matrices_sparse(...) Produces sparse boundary matrices for Betti number computations.

approximate_smith_normal_form(matrix) Approximates the Smith Normal Form of a matrix.

calculate_betti_numbers_and_euler_characteristic(...) Calculates Betti numbers and Euler characteristics.

generate_random_binary_matrix(size, num_black_pixels) Creates a binary matrix with randomly distributed black pixels.

Performance Insights
The computational performance is directly related to the number of black pixels in the matrix. This repository includes a detailed log-log analysis to visualize the relationship.

Contributing
Contributions are welcome! Feel free to fork the repository, make changes, and submit a pull request. For major changes, please open an issue to discuss your proposed improvements.

License
This project is licensed under the MIT License. See the LICENSE file for details.

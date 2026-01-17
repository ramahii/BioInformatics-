import numpy as np
from typing import Union, List, Tuple

class PredictionModel:
    """
    A class to perform discrete step predictions using a transition matrix
    and an initial vector (commonly used in Markov chains, population dynamics, etc.)
    """
    
    def __init__(self, matrix: np.ndarray, initial_vector: np.ndarray):
        """
        Initialize the prediction model.
        
        Args:
            matrix: A square transition matrix (n x n)
            initial_vector: Initial state vector (n,)
        """
        if matrix.shape[0] != matrix.shape[1]:
            raise ValueError("Matrix must be square")
        
        if matrix.shape[0] != len(initial_vector):
            raise ValueError("Matrix and vector dimensions must match")
        
        self.matrix = matrix
        self.initial_vector = initial_vector
        self.states = [initial_vector.copy()]
    
    def predict(self, steps: int = 5) -> List[np.ndarray]:
        """
        Predict the state at each step.
        
        Args:
            steps: Number of discrete steps to predict (default: 5)
        
        Returns:
            List of state vectors for each step (including the initial state)
        """
        current_state = self.initial_vector.copy()
        
        for _ in range(steps):
            # Multiply matrix by current state to get next state
            current_state = self.matrix @ current_state
            self.states.append(current_state.copy())
        
        return self.states
    
    def print_results(self) -> None:
        """Print the prediction results in a formatted way."""
        print("=" * 70)
        print("DISCRETE STEP PREDICTION RESULTS")
        print("=" * 70)
        print(f"\nInitial Vector (Step 0):")
        print(self.states[0])
        print(f"\nMatrix:\n{self.matrix}\n")
        
        for step in range(1, len(self.states)):
            print(f"Step {step}:")
            print(self.states[step])
            print()


def main():
    # Example 1: Simple 2x2 matrix (Markov chain example)
    print("\n" + "=" * 70)
    print("EXAMPLE 1: 2x2 Transition Matrix")
    print("=" * 70)
    
    matrix_1 = np.array([
        [0.7, 0.3],
        [0.4, 0.6]
    ])
    
    vector_1 = np.array([1.0, 0.0])
    
    model_1 = PredictionModel(matrix_1, vector_1)
    model_1.predict(5)
    model_1.print_results()
    
    
    # Example 2: 3x3 matrix (Population dynamics example)
    print("\n" + "=" * 70)
    print("EXAMPLE 2: 3x3 Transition Matrix")
    print("=" * 70)
    
    matrix_2 = np.array([
        [0.5, 0.2, 0.1],
        [0.3, 0.6, 0.2],
        [0.2, 0.2, 0.7]
    ])
    
    vector_2 = np.array([100, 50, 30])
    
    model_2 = PredictionModel(matrix_2, vector_2)
    model_2.predict(5)
    model_2.print_results()
    
    
    # Example 3: User input option
    print("\n" + "=" * 70)
    print("EXAMPLE 3: Custom Input")
    print("=" * 70)
    
    try:
        size = int(input("Enter matrix size (n for nxn matrix): "))
        
        print(f"\nEnter the {size}x{size} matrix row by row:")
        print("(Enter values separated by spaces)")
        
        matrix_3 = []
        for i in range(size):
            row = list(map(float, input(f"Row {i+1}: ").split()))
            if len(row) != size:
                raise ValueError(f"Row must have {size} elements")
            matrix_3.append(row)
        
        matrix_3 = np.array(matrix_3)
        
        print(f"\nEnter the initial vector ({size} values separated by spaces):")
        vector_3 = np.array(list(map(float, input("Vector: ").split())))
        
        if len(vector_3) != size:
            raise ValueError(f"Vector must have {size} elements")
        
        model_3 = PredictionModel(matrix_3, vector_3)
        model_3.predict(5)
        model_3.print_results()
        
    except ValueError as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
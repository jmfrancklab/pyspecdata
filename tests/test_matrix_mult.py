import unittest
from pyspecdata import *
import numpy as np

class TestMatrixMultiplication(unittest.TestCase):
    def test_matrix_multiplication(self):
        # Create two matrices
        matrix1 = nddata(np.array([1, 2, 3, 4]), (2, 2), ['x', 'y'])
        matrix2 = nddata(np.array([5, 6, 7, 8]), (2, 2), ['y', 'z'])

        # Perform matrix multiplication
        result = matrix1 @ matrix2

        # Expected result
        expected_result = nddata(matrix1.data @ matrix2.data, (2, 2), ['x', 'z'])

        # Assert that the result is as expected
        np.testing.assert_array_almost_equal(result.data, expected_result.data)
if __name__ == '__main__':
    unittest.main()

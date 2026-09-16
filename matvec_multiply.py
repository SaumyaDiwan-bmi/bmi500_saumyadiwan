import random

# create a function to compute the dot product of two vectors using a for loop
# LLM Model Used: Claude Sonnet 5 with Medium thinking
# Prompt used: i have a file matvec_test.py. i need to generate code based on the comments
# here is the first comment  import random
# create a function to compute the dot product of two vectors using a for loop
# add comments for the selected function

def dot_product(vec1, vec2):
    """
    Compute the dot product of two vectors using a for loop.
    
    Args:
        vec1 (list): First vector (list of numbers)
        vec2 (list): Second vector (list of numbers)
    
    Returns:
        float/int: The dot product of vec1 and vec2
    """
    # vectors must be the same length to compute a dot product
    if len(vec1) != len(vec2):
        raise ValueError("Vectors must be the same length")

    # accumulator for the running sum of element-wise products
    result = 0

    # loop through each index, multiply corresponding elements, and add to result
    for i in range(len(vec1)):
        result += vec1[i] * vec2[i]

    return result

# create a function to compute the matrix-vector product using the dot_product function
# LLM Model Used: Claude Sonnet 5 with Medium thinking
# Prompt used: here is the second comment # create a function to compute the matrix-vector product using the dot_product function
# add comments for the selected function
def matvec_product(matrix, vec):
    """
    Compute the matrix-vector product using the dot_product function.
    
    Args:
        matrix (list of lists): The matrix, where each inner list is a row
        vec (list): The vector to multiply by
    
    Returns:
        list: The resulting vector from the matrix-vector multiplication
    """
    # number of columns in the matrix must match the length of the vector
    if len(matrix[0]) != len(vec):
        raise ValueError("Number of columns in matrix must match vector length")

    # accumulator list to hold the result of each row's dot product
    result = []

    # loop through each row in the matrix
    for row in matrix:
        # compute the dot product of the current row with the vector
        row_result = dot_product(row, vec)
        # append the result of this row to the output vector
        result.append(row_result)

    return result


# create a main function to test the matrix-vector product function using randomly generated data of size 1000x1000
# LLM Model Used: Claude Sonnet 5 with Medium thinking
# Prompt used: Here is the third comment # create a main function to test the matrix-vector product function using randomly generated data of size 1000x1000

def main():
    # define the size of the matrix and vector
    size = 1000

    # generate a random 1000x1000 matrix with values between 0 and 1
    matrix = [[random.random() for _ in range(size)] for _ in range(size)]

    # generate a random vector of length 1000 with values between 0 and 1
    vec = [random.random() for _ in range(size)]

    # compute the matrix-vector product using our function
    result = matvec_product(matrix, vec)

    # print a small sample of the result to confirm it ran correctly
    print("First 5 elements of the result:", result[:5])
    print("Length of result vector:", len(result))


# run the main function when this script is executed directly
if __name__ == "__main__":
    main()

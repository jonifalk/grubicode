import random
import time

__author__ = "8819989, Falk"
__version__ = "1.0.0"


def generate_random_strings(n_strings: int, len_strings: int) -> list:
    """Generates a list of random DNA sequences."""
    chars = ["A", "T", "C", "G"]
    return [''.join(random.choices(chars, k=len_strings))
            for _ in range(n_strings)]


def benchmark(name: str, func, *args, iterations: int = 1000):
    """Utility function to measure execution time using a loop."""
    start_time = time.perf_counter()
    for _ in range(iterations):
        func(*args)
    end_time = time.perf_counter()
    elapsed = end_time - start_time
    print(
        f"{name: <25} | Iterations: {iterations: <6} "
        f"| Total Time: {elapsed:.5f}s "
        f"| Avg Time: {(elapsed / iterations) * 1e6:.2f} µs")


def p1_version1_loop(sequences):
    """
    Concatenates a list of strings sequentially
    using a for loop and the += operator.
    """
    result = ""
    for s in sequences:
        result += s
    return result


def p1_version2_join(sequences):
    """
    Concatenates a list of strings using the optimized string .join() method.
    """
    return "".join(sequences)


def run_problem_1_experiments():
    """
    Runs benchmarking experiments to compare list concatenation methods.

    Compares the performance of using a "+=" loop versus the ".join()" method
    across two scenarios: many short strings and fewer long strings.
    Returns:
        None: Executes benchmarks and outputs results to the console.
    """
    print("Problem 1: List Concatenation")

    # Scenario A: Many short strings
    seqs_many_short = generate_random_strings(1000, 10)
    print("\nScenario 1A: 1000 strings of length 10")
    benchmark("Version 1 (+= loop)", p1_version1_loop,
              seqs_many_short, iterations=1000)
    benchmark("Version 2 (join)", p1_version2_join,
              seqs_many_short, iterations=1000)

    # Scenario B: Fewer long strings
    seqs_few_long = generate_random_strings(10, 10000)
    print("\nScenario 1B: 10 strings of length 10000")
    benchmark("Version 1 (+= loop)", p1_version1_loop,
              seqs_few_long, iterations=1000)
    benchmark("Version 2 (join)", p1_version2_join,
              seqs_few_long, iterations=1000)
    print("\n")


def p2_version1_add(seqs):
    """Concatenates exactly two strings using the + operator."""
    return seqs[0] + seqs[1]


def p2_version2_format(seqs):
    """Concatenates exactly two strings using the string .format() method."""
    return "{}{}".format(seqs[0], seqs[1])


def p2_version3_fstring(seqs):
    """Concatenates exactly two strings using Python f-strings."""
    return f"{seqs[0]}{seqs[1]}"


def run_problem_2_experiments():
    """
    Runs benchmarking experiments to compare two-string concatenation methods.

    Compares the performance of the "+" operator, ".format()", and f-strings
    across two scenarios: concatenating short strings versus very long strings.

    Returns:
        None: Executes benchmarks and outputs results to the console.
    """
    print("Problem 2: Two-String Concatenation")

    # Scenario A: Short strings
    seqs_short = generate_random_strings(2, 25)
    print("\nScenario 2A: 2 strings of length 25")
    benchmark("Version 1 (+ operator)", p2_version1_add,
              seqs_short, iterations=100000)
    benchmark("Version 2 (.format)", p2_version2_format,
              seqs_short, iterations=100000)
    benchmark("Version 3 (f-string)", p2_version3_fstring,
              seqs_short, iterations=100000)

    # Scenario B: Long strings
    seqs_long = generate_random_strings(2, 50000)
    print("\nScenario 2B: 2 strings of length 50000")
    benchmark("Version 1 (+ operator)", p2_version1_add,
              seqs_long, iterations=100000)
    benchmark("Version 2 (.format)", p2_version2_format,
              seqs_long, iterations=100000)
    benchmark("Version 3 (f-string)", p2_version3_fstring,
              seqs_long, iterations=100000)


if __name__ == "__main__":
    run_problem_1_experiments()
    run_problem_2_experiments()

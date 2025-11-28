## Smith–Waterman Anti-Diagonal Optimized Implementation

This project provides a high-performance implementation of the Smith–Waterman local sequence alignment algorithm using an anti-diagonal dynamic programming (DP) strategy. The code focuses on low-level performance optimization through cache-friendly memory layouts, 64-byte aligned buffers, branch-aware inner loops, and optional OpenMP parallelization for large sequence sizes.

The implementation achieves significant speedups compared to naive DP approaches and reports performance in GFLOPS for benchmarking.

### Features

- Anti-diagonal wavefront traversal for improved data locality

- 64-byte aligned memory allocation for SIMD-friendly access patterns

- Efficient buffer rotation to avoid large DP matrix allocations

- OpenMP parallelization triggered automatically for large inputs

- Configurable scoring (match, mismatch, gap)

- Built-in random DNA sequence generator for testing

- Performance benchmarking with wall-clock timing and GFLOPS estimation


### Build Instructions
Compile with GCC
```
make
```
this will compile the c file

### Usage

Run the program with a single argument specifying the sequence length  (default 1000) and number of threads (default=8):
```
./run.sh baseline
./run.sh optimized
```
baseline will run the python baseline model, while optimized will run the optimized code
Example output:
```
Running baseline python implementation
N=1000 P=8 time=0.449361 seconds
```
```
Using OpenMP with 6 threads
N=1000  score=450  time=0.002625 s
Performance: 12.40 GFLOPS
```
### How It Works

The algorithm processes the DP matrix along anti-diagonals (i + j = constant), allowing each diagonal to be computed independently using only the previous two diagonals. Three aligned buffers are reused and rotated each iteration to avoid storing the entire matrix.

Automatic OpenMP parallelization is applied for large diagonals to distribute work across threads.



These can be adjusted at compile time or edited directly in the source file.

### License

MIT License.

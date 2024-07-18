# uf_decoder

A modified union-find decoder for CSS codes. The decoder makes a simplification as it avoids tracking cluster boundaries and uses instead breadth-first graph traversal of the Tanner graph of the code. The obtained thresholds are on par with the original union-find decoder from https://arxiv.org/abs/1709.06218. Note that this repository is a preliminary version of the decoder and breaking changes may occur in future versions.

## features
- Fast decoding of topological codes with periodic boundary conditions (Algorithms 1, 2, 4)
- Decoding of general qLDPC codes (Algorithm 3)
- Include qubit erasure/loss into the decoding

## usage

The ```C```-code can either be used via one of the examples in the folder ```examples``` or interfaced from ```Python```. For the latter option, you can find examples in the folder ```py_wrapper```. Before usage, the code needs to be compiled with the following commands:

```
$ meson setup build --buildtype=release
$ cd build/
$ ninja
```
For interfacing with ```Python```, you will use the library ```libSpeedDecoder.so``` via ```run_from_py.py```. In our implementation, data qubits and syndromes share the same data structure with syndromes being at the beginning of the array, followed by data qubits. We may change this in a later version to simplify the interfacing from ```Python``` and improve the time and space complexity.

## source code
The implementation of Algorithms 1, 2, 4 can be found in ```src/decoder_main.c```. The implementation of Algorithm 3 is in ```src/decoder_ldpc.c``` with ```src/stabilizer_main.c``` containing the used Gaussian elimination with bitwise-logic.

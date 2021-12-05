About the examples
==================
The examples have been designed with the following issues in mind:
- they run from the Julia REPL
- each example is a Julia module named similar to the basename of the example file.
- an example can be used as the starting point for a project 
- the examples at the same time comprise the test suite for SPHLib.


## Running the examples
In order to run `ExampleXXX`, peform the following steps:

- Download the example file (e.g. via the source code link at the top)
- Call Julia with  an Julia environment which contains `SPHLib.jl`
- `include("ExampleXXX.jl")`
- Run the example via `ExampleXXX.main()`

Due to the encapsulation into modules, you can load as many examples as you like.


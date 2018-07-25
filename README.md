## Basil ##

Basil (for BASIs List) enumerates bases of hyperplane arrangements up to 
symmetries. It depends on David Avis' LRS library to represent the hyperplane 
arrangements and Thomas Rehn's permlib for the permutation calculations.

### Building ###

Run `make basil`

#### Dependencies ####

* Boost (including program_options)
* GMP
* doxygen (for documentation generation)

### Documentation ###

#### Program Options ####

```
  --arrangement-pivot      Makes Basil pivot as if the input was an
                           arrangement.
  --assume-no-symmetry     Forces Basil to assume there is no symmetry in the
                           input.
  --generate-symmetry      Forces Basil to generate a new symmetry group
  --show-all-dicts         Show all intermediate dictionaries in the search
                           tree.
  --no-fixed-plane         Do not fix the x_0 = 1 plane for automorphism
                           calculations. This may result in spurious
                           automorphisms for certain instances (notably those
                           which have rows where x_0 = 0).
  --fund-domain-lim arg    Maximum number of constraints to include in the
                           fundamental domain [default 0]
  --gram arg               Gram matrix generation to use: 'begin' to use the
                           gram matrix from the input file, 'Q' to use the
                           Q-matrix metric for Gram hashing [default],
                           'no-augment' to use the Q-matrix metric without
                           row-augmenting the input first (only works for input
                           matrices of full rank), 'Euclidean' to use the
                           Euclidean metric for Gram matrix generation, or
                           'no-norm' to use the Euclidean metric without
                           normalizing the row vectors of the matrix to the
                           same norm (saves expensive normalization
                           calculations, at the possible expense of not finding
                           all symmetries)
  --no-gram                Deactivates Gram matrix hashing.
  --debug-gram             Print gram vectors for vertices/rays/cobases that
                           are printed
  --stab-search            Activate cobasis stabilizer search (not
                           reccommended, stabilizer computation costs more than
                           it saves)
  --print-basis arg        Print the number of cobases found and running time
                           every n cobases.
  --print-new              Print the added {cobasis,vertex,ray} when printing a
                           status message
  --print-ray arg          Print the number of cobases found and running time
                           every n cobases.
  --print-vertex arg       Print the number of cobases found and running time
                           every n cobases.
  --print-each arg         Convenience for print-{basis,ray,vertex} with the
                           given parameter. If any of the others are given,
                           they take precedence
  --print-trace            Print the full trace of the DFS (warning: very
                           verbose).
  -p [ --preprocess ]      Do not DFS, simply do preprocessing work, and print
                           normalized input file to output stream
  -v [ --verbose ]         Shorthand for --print-interval=1, --print-new. Those
                           settings, if supplied, will take precedence.
  -i [ --input-file ] arg  Input file name. Standard input if none supplied;
                           may also be supplied as first positional argument.
  -m [ --matrix-file ] arg Matrix file name. Alias for --input-file.
  -g [ --group-file ] arg  File to read symmetry group from - overrides any
                           supplied in the input file.
  -o [ --output-file ] arg Output file name. Standard output if none supplied;
                           may also be supplied as second positional argument.
  -h [ --help ]            Produce help message.
```

#### Input Format ####

(based on LRS)

```
[name]
[(H|V|A)-representation]
[linearity <k> < k linearity indices >]
[< other lrs options >]
begin
<n> <d> rational
< n * d whitespace-delimited data values >
end
[symmetry begin
 {<comma-delimeted cycles of whitespace delimeted elements>}
 symmetry end]
[gram auto
|gram Q
|gram augmented
|gram Euclid
|gram inexact
|gram begin
 < n * n whitespace-delimited integers >
 gram end]
[< other lrs options >]
```

#### Class Documentation ####

Run `make doc` (requires doxygen)
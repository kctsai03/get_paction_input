## Usage Instructions

get_paction_inputs take the following parameters and returns a csv with genotype proportions
- decifer output file
- decifer input file (copy number clone proportions)
- the mutation 
- sample_to_pt (i mention what this is in the comments in the code)

get_mut_tree takes the following parameters and returns a csv with the genotype tree edges
- decifer output file
- the mutation
```
An example of using this function is 

test = get_paction_inputs(k14_output, best_input,'10.52573509.G.GA', sample_to_pt)
x = get_mut_tree(k14_output, '10.52573509.G.GA')
```
Output example

test_mut_1.csv

| genotypes                 |s0                | s1        |
| ----------                |----------------  | --------- |
|1\|1\|0         |0.76945265                   |0.362147972764937|
|2\|2\|0	       | 0.2063730203047980        |0.324516934158279|
| 2\|2\|1           |0.0241743346881351        |0.238998250656163|
| 2\|1\|0          |0                          |0.0743368424206201|


test_mut_1_tree.csv

|edge_start| edge_end|
|--------|-----|
| 1\|1\|0               | 2\|2\|0 | 
| 2\|2\|0  |2\|2\|1           |
| 2\|2\|1|2\|1\|0            | 


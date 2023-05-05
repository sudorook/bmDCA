# Example

There are instructions for using bmDCA to train a model on a Kunitz domain
alignment (PF00014, from Pfam database).

## Step 1: Process the input MSA with `process_msa`

Take the provided raw FASTA-formatted alignment and process it to produce a
numerical alignment. The following command will remove sequences with greater
than 20% gaps (`-g`), positions with over 20% gaps (`-G`), and down-weight
sequences with greater than 80% sequence identity (`-t`).

```sh
process_msa -i PF00014_full.fasta -g 0.2 -G 0.2 -t 0.8 -o PF00014_numerical.txt -O PF00014_weights.txt
```

The numerical format alignment is saved as `PF00014_numerical.txt` and the
corresponding sequence weight for each sequence is in `PF00014_weights.txt`.

## Step 2: Run `bmdca`

Use the new numerical alignment and sequence weights to train the model. Run:

```sh
bmdca -n PF00014_numerical.txt -w PF00014_weights.txt -c bmdca.conf -d results
```

The provided `bmdca.conf` file will run for 500 iterations and use the
stochastic gradient descent with momentum to learn the parameters to the Potts
model. Outputs will be saved in the `results/` directory.

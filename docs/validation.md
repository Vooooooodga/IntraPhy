# Software testing and scientific interpretation

Run the unit and raw-input integration suite with:

```bash
python -m unittest discover -s tests -v
```

The suite covers software contracts, numerical behavior, raw FASTA/GFF input
handling, and integration with MAFFT and minimap2. Those external executables
are required for the standard raw-input integration cases. Adapter interface
tests may mock external commands; these tests do not assess the behavior of an
installed external program.

Passing software tests does not establish biological accuracy or statistical
calibration. Inference remains conditional on the supplied annotations, mapped
catalogue, observation rules, and species tree. Empirical accuracy, error rates,
and finite-sample calibration require independent truth or study designs.

## `MeSS` applications

You might want to use `MeSS` when :

- Testing a sequencing run's detection limit:

Let's say, for example, you have a stahpylococcus aureus that is present at roughly 0.01% abundance, and you want to know at which sequencing depth you can detect it. With `MeSS` you can easily make samples with mixes of genomes at different abundances and sequencing depths for a selected sequencing technology.

- Benchmarking taxonomic classifiers:

`MeSS` includes read simulators that generate gold standard bam files. This means that you know the origin of every simulated read, and thus know the true composition of the sample. This allows comparing the results of tools like [centrifuge](https://github.com/DaehwanKimLab/centrifuge) or [kraken2](https://github.com/DerrickWood/kraken2) against a ground truth dataset.

- Testing the impact of sequencing run parameters on genome assembly

Let's suppose that you have a bacterial isolate that you plan to sequence. However, you are not sure about which sequencing technology to use, at which depth, etc ...
You can use `MeSS` to easily to generate reads at different depth, with different sequencing technologies and error profiles (ONT 10.4.1, PacBio hifi, illumina miseq...). Simulated reads can be assembled, and you an evaluate the impact of sequencing variables on your assembly.

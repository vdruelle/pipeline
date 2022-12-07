import click
import os
import numpy as np
from Bio import Phylo, SeqIO


@click.group()
def cli():
    pass

@cli.command()
@click.argument("newick_tree", type=click.Path(exists=True))
@click.argument("clusters_folder", type=click.Path(exists=True))
@click.argument("output_name")
def newick_to_txt(newick_tree, clusters_folder, output_name):
    """
    Converts newick tree to txt for viewing.
    """
    tree = Phylo.read(newick_tree, "newick")
    tips = tree.get_terminals()
    tree_length = max([tree.distance(tip) for tip in tips])
    if tree_length < 0.001:
        with open(clusters_folder + "/contigs_tree.txt", "w") as f:
                f.write("Tree too small for plotting, which means the contig is fine.")
    else:
        os.system(f"nw_display {newick_tree} > {output_name}")



@cli.command()
@click.argument("clusters_folder", type=click.Path(exists=True))
@click.argument("tree_file", type=click.Path(exists=True))
@click.argument("output", type=click.Path(exists=False))
@click.argument("max_length", type=float)
def check_contigs(clusters_folder, tree_file, output, max_length):
    """
    Check generated contigs to decide whether the pipeline should continue automatically or if it requires
    human intervention.
    """
    ### Check number of contigs created
    folder_count = 0
    for files in os.listdir(clusters_folder):  # loop over all files
        if os.path.isdir(os.path.join(clusters_folder, files)):  # if it's a directory
            folder_count += 1
    
    print(f"Numbers of clusters created: {folder_count}")
    if folder_count == 0:
        print("No cluster detected, cannot continue.")
        with open(output, "w") as f:
                f.write("No cluster detected, cannot continue. Investigate the error.")
    elif folder_count > 1:
        print("More than one cluster detected, manual intervention needed before next part of pipeline.")
        with open(output, "w") as f:
                f.write("More than one cluster detected, remove the cluster you do not wish to reconstruct.")
    elif folder_count == 1:
        print("Only one cluster, continuing with inspection of contigs tree.")
        
        # Getting tree length
        tree = Phylo.read(tree_file, "newick")
        tips = tree.get_terminals()
        tree_length = max([tree.distance(tip) for tip in tips])
        if tree_length > max_length:
            print("Assemblies in the contig are too diverged, manual intervention needed before next part of \
                   pipeline.")
            with open(output, "w") as f:
                f.write("Assemblies in the contig are too diverged, manual intervention might be needed \
                         before next part of the pipeline.")
                f.write("Look at the contig tree and remove the contigs you don't want to include.")
                f.write("Then proceed with the second part of the pipeline.")
        else:
            print("Assemblies are similar enough, you can run the second part of the pipeline.")
            with open(output, "w") as f:
                f.write("All good, the next par of the pipeline can be run.")

@cli.command()
@click.argument("combined", type=click.Path(exists=True))
def rotate(combined):
    """
    Rotates the sequence so that the starting point are the same.
    """
    records = list(SeqIO.parse(combined,"fasta"))
    illumina = records[0]
    nanopore = records[1]
    nanopore_complement = records[1].reverse_complement()

    klen = 30
    found_start = False
    while not found_start:
        query = illumina.seq[:klen]
        tmp = nanopore.seq.find(query)
        tmp2 = nanopore_complement.seq.find(query)
        if tmp != -1 or tmp2 != -1:
            if tmp != -1:
                nanopore.seq = nanopore.seq[tmp:] + nanopore.seq[:tmp]
            elif tmp2 != -1:
                print("Nanopore sequence had to be reverse-complemented.")
                nanopore.seq = nanopore_complement.seq[tmp2:] + nanopore_complement.seq[:tmp2]
            found_start = True
            print(f"Nanopore sequence was rotated succesfully with klen={klen}.")
        else:
            if klen > 10:
                klen -= 1
            else:
                break
    
    if not found_start:
        print("Nanopore sequence was not rotated succesfully.")

    with open(combined, "w") as f:
        SeqIO.write([illumina, nanopore], f, "fasta")


@cli.command()
@click.argument("alignment", type=click.Path(exists=True))
def statistics(alignment):
    """
    Compute quick statistics on the accuracy of the nanopore assembly.
    """
    records = list(SeqIO.parse(alignment,"fasta"))
    illumina = np.array(records[0])
    nanopore = np.array(records[1])

    diff_pos = np.where(illumina != nanopore)[0]
    nb_diff = diff_pos.shape[0]

    print(f"Sequences differ at {nb_diff} nucleotide positions.")
    print(f"These posistions are:")
    print(diff_pos)


if __name__ == '__main__':
    cli()

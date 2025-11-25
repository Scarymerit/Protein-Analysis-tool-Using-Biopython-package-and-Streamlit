# ===============================
# Install required packages
# ===============================


# ===============================
# Imports
# ===============================
import os, tempfile, random, requests
from Bio import SeqIO, motifs, Phylo, pairwise2, Entrez, ExPASy
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBParser, PPBuilder
from IPython.display import display, Image
import nglview as nv
import matplotlib.pyplot as plt
import ipywidgets as widgets

# Create results folder
os.makedirs("results", exist_ok=True)


# ===============================
# 1. Upload PDB file
# ===============================
def upload_pdb():
    from google.colab import files
    print("\n=== Upload PDB File ===")
    uploaded = files.upload()
    pdb_file = list(uploaded.keys())[0]
    print(f"PDB file uploaded: {pdb_file}")
    return pdb_file


# ===============================
# 2. Extract peptide sequences from PDB
# ===============================
def extract_sequences_from_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)
    ppb = PPBuilder()
    sequences = []
    for model in structure:
        for chain in model:
            peptides = ppb.build_peptides(chain)
            for i, pep in enumerate(peptides):
                sequences.append(SeqRecord(Seq(str(pep.get_sequence())), id=f"{chain.id}_{i}"))
    print(f"Extracted {len(sequences)} peptide sequences from PDB")
    return sequences


# ===============================
# 3. Reference-based MSA
# ===============================
def manual_reference_msa(sequences):
    if not sequences: return None
    reference = max(sequences, key=lambda s: len(s.seq))
    ref_len = len(reference.seq)
    aligned_records = [SeqRecord(reference.seq, id=reference.id)]
    for seq in sequences:
        if seq.id == reference.id: continue
        alignment = pairwise2.align.globalms(reference.seq, seq.seq, 2, -1, -0.5, -0.1)[0]
        aligned_seq = alignment[1]
        if len(aligned_seq) < ref_len:
            aligned_seq += "-" * (ref_len - len(aligned_seq))
        elif len(aligned_seq) > ref_len:
            aligned_seq = aligned_seq[:ref_len]
        aligned_records.append(SeqRecord(Seq(aligned_seq), id=seq.id))
    msa = MultipleSeqAlignment(aligned_records)
    print(f"MSA completed: {len(msa)} sequences, length {msa.get_alignment_length()} bp")
    return msa


# ===============================
# 4. Motif Discovery with conserved positions
# ===============================
def motif_discovery(msa):
    if not msa: return [], [], [], {}
    seq_objs = [record.seq for record in msa]
    motif_obj = motifs.create(seq_objs)
    logo_file = "results/motif_logo.png"
    motif_obj.weblogo(logo_file)
    print(f"Motif logo saved as {logo_file}")

    # Conserved positions
    alignment_length = msa.get_alignment_length()
    conserved_pos = [i for i in range(alignment_length) if len(set(msa[:, i].replace("-", ""))) == 1]

    # Per-chain motif mapping
    motif_mapping = {}
    for record in msa:
        motif_mapping[record.id] = [i for i in range(alignment_length) if i in conserved_pos]

    return [motif_obj], [logo_file], conserved_pos, motif_mapping


# ===============================
# 5. Phylogenetic Trees
# ===============================
def phylogenetic_tree(msa):
    if not msa or len(msa) < 2:
        print("Not enough sequences for phylogenetic tree (needs >=2). Skipping tree construction.")
        return None, None, None

    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(msa)
    constructor = DistanceTreeConstructor()
    upgma_tree = constructor.upgma(dm)
    nj_tree = constructor.nj(dm)

    # Save trees
    Phylo.write(upgma_tree, "results/upgma_tree.newick", "newick")
    Phylo.write(nj_tree, "results/nj_tree.newick", "newick")
    for tree, fname in zip([upgma_tree, nj_tree], ["results/upgma_tree.png", "results/nj_tree.png"]):
        plt.figure(figsize=(12, 6))
        Phylo.draw(tree, do_show=False)
        plt.savefig(fname)
        plt.close()
    print("Phylogenetic trees saved as Newick + PNG")
    return upgma_tree, nj_tree, dm


# ===============================
# 6. Bootstrap Phylogenetic Trees
# ===============================
def bootstrap_phylo(msa, n_bootstrap=100):
    if not msa or len(msa) < 2:
        print("Not enough sequences for bootstrap. Skipping bootstrap.")
        return []
    calc = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor()
    trees = []
    for _ in range(n_bootstrap):
        sampled_cols = [random.randint(0, msa.get_alignment_length() - 1) for _ in range(msa.get_alignment_length())]
        new_records = []
        for rec in msa:
            new_seq = ''.join([rec.seq[i] for i in sampled_cols])
            new_records.append(SeqRecord(Seq(new_seq), id=rec.id))
        new_msa = MultipleSeqAlignment(new_records)
        dm = calc.get_distance(new_msa)
        tree = constructor.nj(dm)
        trees.append(tree)
    print(f"Bootstrap generated {len(trees)} trees")
    return trees


# ===============================
# 7. Interactive 3D Motif Visualization
# ===============================
def interactive_3d_motifs(pdb_file, motif_mapping):
    view = nv.show_file(pdb_file)
    view.add_cartoon()
    if not motif_mapping:
        display(view)
        return view

    chain_selector = widgets.Dropdown(
        options=[chain for chain in motif_mapping.keys()],
        description='Select Chain:',
        value=list(motif_mapping.keys())[0]
    )
    out = widgets.Output()

    def highlight_chain(change):
        out.clear_output()
        view.clear_representations()
        view.add_cartoon()
        selected_chain = change['new']
        positions = motif_mapping[selected_chain]
        if positions:
            selection = f":{selected_chain} and ({' or '.join([str(p + 1) for p in positions])})"
            view.add_cartoon(selection=selection, color='red')
        view.center()
        with out: display(view)

    chain_selector.observe(highlight_chain, names='value')
    display(chain_selector, out)
    highlight_chain({'new': chain_selector.value})
    return view


# ===============================
# 8. Hydrophobicity Mapping
# ===============================
def hydrophobicity_3d(pdb_file, pdb_sequences):
    if not pdb_sequences: return
    hydrophobic_residues = "AILMFWV"
    view = nv.show_file(pdb_file)
    view.add_cartoon()
    for record in pdb_sequences:
        chain_id = record.id
        selection = []
        for i, aa in enumerate(record.seq):
            if aa.upper() in hydrophobic_residues:
                selection.append(str(i + 1))
        if selection:
            view.add_cartoon(selection=f":{chain_id} and ({' or '.join(selection)})", color='orange')
    view.center()
    display(view)
    return view


# ===============================
# 9. Motif Frequency Plot
# ===============================
def motif_frequency_plot(msa):
    if not msa: return
    alignment_length = msa.get_alignment_length()
    freq = [len(set(msa[:, i].replace("-", ""))) for i in range(alignment_length)]
    plt.figure(figsize=(12, 4))
    plt.bar(range(1, alignment_length + 1), freq)
    plt.xlabel("Position")
    plt.ylabel("Unique residues (motif frequency)")
    plt.title("Motif Frequency Across Alignment")
    plt.show()


# ===============================
# 10. Fetch GenBank/SwissProt
# ===============================
def fetch_genbank_swissprot(ids):
    if not ids: return [], []
    Entrez.email = "your.email@example.com"
    gb_records, sp_records = [], []
    for id in ids:
        try:
            handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
            gb_records.append(SeqIO.read(handle, "genbank"))
            handle.close()
        except:
            print(f"Warning: Could not fetch GenBank ID {id}")
        try:
            handle = ExPASy.get_sprot_raw(id)
            sp_records.append(SeqIO.read(handle, "swiss"))
            handle.close()
        except:
            print(f"Warning: Could not fetch SwissProt ID {id}")
    return gb_records, sp_records


# ===============================
# 12. Visual Summary
# ===============================
def visual_summary(msa, motif_images, conserved_pos, motif_mapping, upgma_tree, nj_tree, pdb_file, pdb_sequences,
                   gb_records=[], sp_records=[]):
    print("\n=== VISUAL SUMMARY ===")
    # MSA stats
    print("\nMSA Stats:")
    for r in msa:
        seq_str = str(r.seq)
        gaps = seq_str.count('-')
        print(f"{r.id} | Length: {len(seq_str)} | Gaps: {gaps} | %Gaps: {gaps / len(seq_str) * 100:.2f}")
    # Motif logos
    print("\nMotif Logos:")
    for img in motif_images: display(Image(img))
    print(f"Conserved positions: {conserved_pos}")
    # Motif frequency plot
    motif_frequency_plot(msa)
    # Phylo toggle
    if upgma_tree and nj_tree:
        tree_selector = widgets.ToggleButtons(options=["UPGMA", "NJ"], description="Tree:")
        out = widgets.Output()

        def show_tree(change):
            out.clear_output()
            with out:
                if change['new'] == "UPGMA":
                    Phylo.draw(upgma_tree)
                else:
                    Phylo.draw(nj_tree)

        tree_selector.observe(show_tree, names='value')
        display(tree_selector, out)
    else:
        print("Phylogenetic tree not available (MSA < 2 sequences).")
    # 3D Motif visualization
    print("\nInteractive 3D Motif Visualization:")
    interactive_3d_motifs(pdb_file, motif_mapping)
    # Hydrophobicity
    print("\n3D Hydrophobic Residue Mapping (orange):")
    hydrophobicity_3d(pdb_file, pdb_sequences)
    # GenBank/SwissProt
    if gb_records:
        print("\nGenBank Records:")
        for gb in gb_records: print(f"ID: {gb.id}, Desc: {gb.description}, Length: {len(gb.seq)}")
    if sp_records:
        print("\nSwissProt Records:")
        for sp in sp_records: print(f"ID: {sp.id}, Desc: {sp.description}, Length: {len(sp.seq)}")


# ===============================
# 13. Main Pipeline
# ===============================
def main_pipeline():
    print("\n" + "=" * 60)
    print("WELCOME TO THE ADVANCED PROTEIN ANALYSIS PIPELINE")
    print("=" * 60)

    # Upload PDB
    pdb_file = upload_pdb()
    # Extract sequences
    pdb_sequences = extract_sequences_from_pdb(pdb_file)
    # MSA
    pdb_msa = manual_reference_msa(pdb_sequences)
    # Motif discovery
    motif_list, motif_images, conserved_pos, motif_mapping = motif_discovery(pdb_msa)
    # Phylogenetic trees
    upgma_tree, nj_tree, dm = phylogenetic_tree(pdb_msa)
    # Optional GenBank/SwissProt
    gb_ids, sp_ids = ["NM_000546"], ["P04637"]
    gb_records, sp_records = fetch_genbank_swissprot(gb_ids + sp_ids)
    # Visual summary
    visual_summary(pdb_msa, motif_images, conserved_pos, motif_mapping, upgma_tree, nj_tree, pdb_file, pdb_sequences,
                   gb_records, sp_records)
    print("\nPIPELINE COMPLETED SUCCESSFULLY!")


# ===============================
# Run Pipeline
# ===============================
if _name == "main_":
    main_pipeline()
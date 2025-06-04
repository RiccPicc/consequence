from Bio import Entrez, SeqIO
import os

Entrez.email = "your_email@example.com"


def build_search_query(gene_name, taxon_name=None, taxid=None):
    """Build a search query for NCBI nucleotide database."""

    parts = [f"{gene_name}[All Fields]"]
    
    if taxid:
        parts.append(f"txid{taxid}[Organism:exp]")
    elif taxon_name:
        parts.append(f"{taxon_name}[Organism]")
    else:
        raise ValueError("Either taxid or taxon_name must be provided.")
    
    query = " AND ".join(parts)
    return query


def search_nucleotide_ids(query, max_results):
    """Search NCBI nucleotide database for IDs matching query."""
    print(f"Searching NCBI with query: {query}")
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]


def download_sequences(ids, out_dir):
    """Download sequences by IDs in GenBank or FASTA format."""
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for seq_id in ids:
        try:
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
            filename = os.path.join(out_dir, f"{seq_id}.fasta")
            with open(filename, "w") as f:
                f.write(handle.read())
            handle.close()
            print(f"Saved: {filename}")
            return handle.read()
        except Exception as e:
            print(f"Error downloading {seq_id}: {e}")


def run_download(gene_name, taxon_name=None, taxid=None, max_results=100, out_dir="."):
    query = build_search_query(gene_name, taxon_name=taxon_name, taxid=taxid)
    ids = search_nucleotide_ids(query, max_results=max_results)
    if not ids:
        print("No sequences found.")
        return None
    else:
        print(f"Found {len(ids)} sequence(s). Downloading...")
        return download_sequences(ids, out_dir=out_dir)
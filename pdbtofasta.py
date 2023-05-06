import re
import os
import shutil
from Bio.PDB import PDBParser, Select, PDBIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import glob
from itertools import chain

aa_dict = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

if os.path.exists("heavy_chains"):
    shutil.rmtree("heavy_chains")

if os.path.exists("light_chains"):
    shutil.rmtree("light_chains")

os.makedirs("heavy_chains")
os.makedirs("light_chains")


class ChainSelector(Select):
    def __init__(self, chain_ids):
        self.chain_ids = chain_ids

    def accept_chain(self, chain):
        return chain.get_id() in self.chain_ids


def extract_chains(pdb_file, chain_ids):
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("temp", pdb_file)
    except ValueError as e:
        print(f"Error processing {pdb_file}: {e}")
        return []

    chains = []

    for chain in structure.get_chains():
        if chain.get_id() in chain_ids:
            seq = ""
            for residue in chain.get_residues():
                if residue.get_resname() in aa_dict:
                    seq += aa_dict[residue.get_resname()]
            matches = re.findall(r"\d+", str(chain.get_full_id()))
            if matches:
                chain_id = matches[0]
            else:
                chain_id = "unknown"
            filename = os.path.basename(pdb_file)
            record_id = f"{filename}_{chain_id}"
            chains.append(SeqRecord(Seq(seq), id=record_id, description=""))

    return chains



pdb_files = glob.glob("all_structures/chothia1/*.pdb")

for pdb_file in pdb_files:
    heavy_chains = extract_chains(pdb_file, ["H"])
    light_chains = extract_chains(pdb_file, ["L"])
    filename = os.path.splitext(os.path.basename(pdb_file))[0]
    SeqIO.write(heavy_chains, f"heavy_chains/{filename}_heavy.fasta", "fasta")
    SeqIO.write(light_chains, f"light_chains/{filename}_light.fasta", "fasta")
import logging
import os

logger = logging.getLogger("whiscy_setup")


def _parse_hssp_proteins(line_buffer):
    """Parses the '## PROTEINS section' from HSSP alignment file"""
    proteins = {}
    for line in line_buffer:
        if line.startswith("  NR.") or line.startswith("##"):
            continue
        # Only get the id and name of the protein in the alignment
        fields = (line[:20]).split(":")
        seq_id = int(fields[0]) - 1
        name = fields[1].strip()[:10]
        proteins[seq_id] = name
    return proteins


def _parse_hssp_alignments(line_buffer, chain_id, num_alignments):
    """Parses the '## ALIGNMENTS section' from HSSP alignment file"""
    alignments = [[] for i in range(num_alignments)]
    # Residues of the chain itself, one per alignment column
    master = []
    first_alignment = 0
    last_alignment = 0
    current_num_alignments = 0
    for line in line_buffer:
        if line.startswith(" SeqNo") or line[12] == "!":
            continue
        if line.startswith("## ALIGNMENTS"):
            fields = (line[13:]).split("-")
            # We are now parsing alignments from first to last specified
            # in the ALINGMENTS header
            first_alignment = int(fields[0]) - 1
            last_alignment = int(fields[1]) - 1
            current_num_alignments = last_alignment - first_alignment + 1
        else:
            if line[12] == chain_id and line[14] != "X":
                # Only collect the chain's residues once, not per alignment block
                if first_alignment == 0:
                    # DSSP writes disulfide-bonded cysteines as lower case letters
                    master.append("C" if line[14].islower() else line[14])
                for i, s in enumerate(line[51 : 51 + current_num_alignments]):
                    # We will convert spaces or dots to -
                    if s == "." or s == " ":
                        s = "-"
                    # We leave residues in minor case as if to not forget insertions
                    alignments[first_alignment + i].append(s)
    alignments = [("".join(s)) for s in alignments]
    return "".join(master), alignments


def hssp_file_to_phylip(hssp_file_name, phylip_file_name, chain_id, master_sequence):
    """Parses an HSSP file and returns a list of the sequences"""
    # The MASTER row is the chain sequence as listed in the HSSP file, so it has
    # one residue per alignment column. It can differ from master_sequence (read
    # from the PDB file), e.g. residues with insertion codes; whiscy_setup maps
    # the PDB residues onto it afterwards.
    # We're only interested in the lenght of the sequence of our given chain_id,
    # SEQLENGHT header gives us the sum of all.
    seqlength = len(master_sequence)
    num_alignments = 0
    num_chains = 0
    parsing = False
    parsing_alignment = False
    line_buffer = []
    parsing_proteins = False
    prot_line_buffer = []
    with open(hssp_file_name, "r") as handle:
        for line in handle:
            line = line.rstrip(os.linesep)
            if line.startswith("NCHAIN"):
                num_chains = int(line.split()[1])
            if line.startswith("NALIGN"):
                num_alignments = int(line.split()[1])

            parsing = seqlength != 0 and num_chains != 0 and num_alignments != 0

            if parsing:
                if line.startswith("## ALIGNMENTS"):
                    parsing_alignment = True

                if line.startswith("## PROTEINS"):
                    parsing_proteins = True

                if line.startswith("##") and "ALIGNMENTS" not in line:
                    parsing_alignment = False

                if line.startswith("##") and "PROTEINS" not in line:
                    parsing_proteins = False

                if parsing_alignment:
                    line_buffer.append(line)

                if parsing_proteins:
                    prot_line_buffer.append(line)

        proteins = _parse_hssp_proteins(prot_line_buffer)
        hssp_master, alignments = _parse_hssp_alignments(
            line_buffer, chain_id.upper(), num_alignments
        )
        seqlength = len(hssp_master)

        all_zero = sum([len(a) for a in alignments]) == 0

        if all_zero:
            raise Exception(
                "Not a single alignment found for chain {}".format(chain_id)
            )

        non_valid = [
            k for k in proteins.keys() if alignments[k].count("-") >= seqlength
        ]
        with open(phylip_file_name, "w") as output_handle:
            # Write header, MASTER also counts
            output_handle.write(
                "{}  {}{}".format(
                    len(proteins) - len(non_valid) + 1, seqlength, os.linesep
                )
            )
            # Write master sequence
            output_handle.write("MASTER    {}{}".format(hssp_master, os.linesep))
            # Write the rest of non null alignments
            for k in sorted(proteins.keys()):
                if k not in non_valid:
                    output_handle.write(
                        "{:10s}{}{}".format(proteins[k], alignments[k], os.linesep)
                    )

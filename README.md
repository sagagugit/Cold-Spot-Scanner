# How to run cold spot scanner 

[Open in Colab](https://colab.research.google.com/github/sagagugit/Cold-Spot-Scanner/blob/main/Cold_Spot_Scanner.ipynb)

About cold-spots in protein-protein interactions

Cold spots are positions in proteins occupied by non-optimal amino acids. At cold-spot positions, several possible mutations lead to affinity enhancement. We identified three possible scenarios how cold-spots could occur:

(I) Cold-spots occur at positions where a WT amino acid does not interact with the binding partner, creating a cavity.

(II) Cold-spots occur at positions where a charged amino acid is buried in hydrophobic environment, exhibiting unfavourable Charge-Hydrophobic (CH) interactions. Eliminating the charged amino acid through a mutation enhances affinity.

(III) Cold-spots occur at positions where two amino acids of the same charge come close together, exibiting unfavorable Same-Charge (SC) interactions. Eliminating one of these charges through a mutation enhances affinity.

Previously, we searched for cold spots in the whole PDB and proved that cold spots are a general feature in evolution of protein-protein iteractions.

About this Colab notebook

This colab notebook allows you to identify cold spots due to the above three scenarios in your favorite protein-protein complex. Once the structure is uploaded, the notebook computes the cold spots, visualizes the cold spots in the structure of the complex, and creates a summary file for download.

Instructions for using this Colab notebook

Two options are possible for uploading the protein complex structure.

1) The complex structure is downloaded directly from the PDB. Please input the "PDB ID" of the Protein complex and chain IDs to define the interface of the protein complex under the sub category "PDB ID and chains". If you mention chains that are not interacting, the program will generate an error message upon execution.

2) The complex structure is uploaded from the user’s computer. To enable users to upload their own complex, kindly remove the comment symbols (#) from all lines in the section labeled "Uploading the complex instead of PDB ID". Once uncommented, the user can upload their desired complex upon execution. Before execution of the program, specify its corresponding four-letter name and chains under the "PDB ID and chains" subcategory. Please ensure that the name assigned to the uploaded complex is limited to a four-letter format, similar to that of a PDB ID.

# Virtual PCR Simulator

🔬 **Virtual PCR Simulator** is an interactive web application built with Streamlit and Biopython that allows users to simulate the Polymerase Chain Reaction (PCR) process **in silico**.  
It accepts DNA sequences via typing or FASTA file upload, accepts forward and reverse primers, predicts all valid PCR products, and visualizes them on a simulated gel electrophoresis image.

---

## 🚀 Features

- Input DNA sequence by typing/pasting or uploading FASTA files (`.fasta`, `.fa`).  
- Specify forward and reverse primers.  
- Detects multiple PCR products based on primer matches and their positions.  
- Displays PCR product details in a scrollable, sortable table:  
  - Product number  
  - Start and end positions  
  - Length (bp)  
  - GC content (%)  
  - Sequence preview (first 50 bp)  
- Simulates gel electrophoresis bands with size annotations.  
- Clean, user-friendly interface with validation and helpful messages.

---

## 🛠️ Installation

pip install:

streamlit
biopython
matplotlib
pandas

RUN THE APP LOCALLY WITH :  streamlit run app.py




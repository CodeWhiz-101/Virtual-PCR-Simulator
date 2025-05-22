import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO
import matplotlib.pyplot as plt
import pandas as pd

st.set_page_config(page_title="Virtual PCR Simulator")
st.title("🔬 Virtual PCR Simulator")
input_method = st.radio("How do you want to provide DNA?", ["Type / Paste it", "Upload FASTA file"])
dna_seq = ""

if input_method == "Type / Paste it":
    dna_seq = st.text_area("🧬 DNA Sequence (5' → 3')", height=150)
else:
    uploaded_file = st.file_uploader("📁 Upload DNA FASTA", type=["fasta", "fa"])
    if uploaded_file is not None:
        try:
            text = uploaded_file.getvalue().decode("utf-8")
            record = SeqIO.read(StringIO(text), "fasta")
            dna_seq = str(record.seq)
            st.success(f"✅ Loaded sequence: {record.id}")
            st.text(f"First 100 bp: {dna_seq[:100]}...")
        except Exception as e:
            st.error(f"❌ Could not read FASTA file: {e}")

fwd_primer = st.text_input("➡️ Forward Primer (5' → 3')")
rev_primer = st.text_input("⬅️ Reverse Primer (5' → 3')")

if st.button("🔁 Run PCR"):
    dna = dna_seq.upper().replace(" ", "").replace("\n", "")
    fwd = fwd_primer.upper()
    revRC = str(Seq(rev_primer).reverse_complement())

    if not dna:
        st.error("❗ Please provide DNA sequence first!")
    elif not fwd or not rev_primer:
        st.error("❗ Please enter both primers!")
    else:

        def find_all(sequence, sub):
            positions = []
            start = 0
            while True:
                start = sequence.find(sub, start)
                if start == -1:
                    break
                positions.append(start)
                start += 1
            return positions

        fwd_positions = find_all(dna, fwd)
        rev_positions = find_all(dna, revRC)
        valid_products = []

        for fwd_pos in fwd_positions:
            for rev_pos in rev_positions:
                if fwd_pos < rev_pos:
                    product = dna[fwd_pos : rev_pos + len(revRC)]
                    valid_products.append((fwd_pos, rev_pos + len(revRC), product))

        if valid_products:
            data = []
            for i, (start, end, product) in enumerate(valid_products, 1):
                gc = 100 * (product.count('G') + product.count('C')) / len(product)
                data.append({
                    "Product #": i,
                    "Start": start,
                    "End": end,
                    "Length (bp)": len(product),
                    "GC Content (%)": f"{gc:.2f}",
                    "Seq (first 50 bp)": product[:50] + ("..." if len(product) > 50 else "")
                })
            df = pd.DataFrame(data)

            st.markdown("### 🧬 Multiple PCR Products")
            st.dataframe(df, height=300)

            st.markdown("### 🧫 Simulated Gel Electrophoresis")
            lengths = [len(p) for (_, _, p) in valid_products]
            fig, ax = plt.subplots(figsize=(2.5, 5))
            max_bp = max(lengths + [1000])
            for i, length in enumerate(lengths):
                ypos = max_bp - length
                ax.plot([0.5], [ypos], marker='s', color='black', markersize=15)
                ax.text(0.6, ypos, f"{length} bp", va='center', fontsize=8)
            ax.set_ylim(0, max_bp)
            ax.set_xlim(0, 1)
            ax.axis('off')
            st.pyplot(fig)

        else:
            st.warning("⚠️ No valid PCR products found with given primers.")


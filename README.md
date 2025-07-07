# 🔬 Molecular Morphing via Beam Search
> A bidirectional molecule transformation pipeline using RDKit, combining core extraction, beam search, and synthetic accessibility optimization.

---

## ✨ Overview

This project performs **molecular morphing** between two molecules given as SMILES strings. It incrementally transforms the **source** molecule into the **target** while:
- Preserving a core scaffold
- Maximizing structural and fingerprint similarity
- Minimizing synthetic complexity

It generates and logs all **valid intermediate structures** into a CSV file with computed metrics.

---

## 🚀 Features

- ✅ Core extraction via **SMARTS**, **Murcko Scaffold**, or **MCS**
- 🔁 Bidirectional morphing from source ⇄ target
- 📈 SA score, Tanimoto similarity, and composite scoring
- 🌱 Growth via beam search and BFS
- 📊 CSV logging of all intermediates
- 🧬 Molecule visualization (grid)

---

## 📦 Installation

```bash
pip install -q rdkit molvs sascorer matplotlib pillow

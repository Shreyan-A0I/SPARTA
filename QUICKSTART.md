# SPARTA Quick Start Guide

## Installation

```bash
# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Running the Application

```bash
streamlit run app.py
```

The app opens at `http://localhost:8501`.

---

## Mode Selection

SPARTA offers two analysis modes:

| Mode | Use Case |
|------|----------|
| **Single MSI** | Compare two *different* metabolites within one tissue section |
| **Comparison** | Compare the *same* metabolite across two MSI datasets (pre/post) |

Switch modes using the radio button at the top of the sidebar.

---

## Single MSI Mode

### Workflow

1. **Load Data**
   - Check "Use local file path" → Enter `data/mouse1.imzML`
   - Click "Load MSI 1"
   - Upload METASPACE CSV (optional but recommended)

2. **Select Metabolites**
   - Choose **Metabolite A** from dropdown
   - Choose **Metabolite B** from dropdown
   - Or enter m/z values manually

3. **Align**
   - Click "🔄 Auto-Align" for automatic centroid alignment
   - Fine-tune with nudge buttons (⬅️⬆️⬇️➡️) if needed

4. **Visualize**
   - **Top Panel**: Side-by-side heatmaps + RGB overlay
   - **Co-loc Map**: Log₂ ratio showing A vs B dominance

5. **Export**
   - Click "📊 Export CSV" for pixel-wise data
   - Click "🖼️ Export PNG" for figure panel

### Interpreting the Co-localization Map

| Color | Meaning |
|-------|---------|
| 🔴 Red | Metabolite A dominant (A >> B) |
| 🔵 Blue | Metabolite B dominant (B >> A) |
| ⚪ White | Balanced / Co-localized (A ≈ B) |
| ⬛ Black | Below threshold (masked) |

---

## Comparison Mode (Pre/Post Analysis)

### Workflow

1. **Switch to Comparison Mode**
   - Select "Comparison Mode" in sidebar

2. **Load MSI 1 (Pre-condition)**
   - Enter path `data/mouse1.imzML` or upload
   - Upload corresponding CSV annotations

3. **Load MSI 2 (Post-condition)**
   - Enter path `data/mouse2.imzML` or upload
   - Upload corresponding CSV annotations

4. **Select Target Metabolite**
   - App shows **common metabolites** found in both CSVs
   - Select the metabolite you want to compare

5. **Align MSIs**
   - Click "🔄 Auto-Align MSIs" to match centroids
   - Adjust with nudge controls if needed

6. **Analyze**
   - **Individual Heatmaps**: Side-by-side Pre vs Post
   - **Ratio Heatmap**: Post/Pre, Pre/Post, or Log₂(Post/Pre)
   - **Overlay**: RGB composite (Red=Pre, Blue=Post)
   - **Difference Map**: Normalized (Post - Pre)
   - **Statistics**: Mean, coverage, fold change
   - **Centroid Mapping**: Spatial shift visualization

7. **Export**
   - CSV with pre/post intensities and ratios
   - PNG panel with all visualizations + histogram

### Ratio Types

| Option | Formula | Best For |
|--------|---------|----------|
| Post/Pre | Post ÷ Pre | Upregulation analysis |
| Pre/Post | Pre ÷ Post | Downregulation analysis |
| Log₂(Post/Pre) | log₂(Post/Pre) | Symmetric fold-change view |

---

## Example Data

The `data/` folder contains example datasets:

```
data/
├── mouse1.imzML        # MSI Dataset 1
├── mouse1.ibd
├── mouse1_pixel_intensities.csv
├── mouse2.imzML        # MSI Dataset 2
├── mouse2.ibd
└── mouse2_pixel_intensities.csv
```

**Quick Test (Single MSI):**
```
Path: data/mouse1.imzML
CSV: data/mouse1_pixel_intensities.csv
```

**Quick Test (Comparison):**
```
MSI 1: data/mouse1.imzML + mouse1_pixel_intensities.csv
MSI 2: data/mouse2.imzML + mouse2_pixel_intensities.csv
```

---

## Testing

```bash
# Run all tests
pytest tests/ -v

# Run specific module
pytest tests/test_centroid_engine.py -v

# Run with coverage report
pytest tests/ --cov=backend --cov-report=html
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| **Import errors** | Run `pip install -r requirements.txt` |
| **File not found** | Ensure `.imzML` and `.ibd` are in same directory |
| **No signal** | Lower SNR floor multiplier (try 1.5–2.0) |
| **Centroid failed** | Check m/z values, try different metabolite |
| **Misaligned images** | Use manual nudge after auto-align |
| **No common metabolites** | Uncheck "Show only common metabolites" |

---

## Keyboard Shortcuts

- The app is fully mouse-driven
- Use Tab to navigate between controls
- Enter to confirm text inputs

---

## Tips

1. **Start with auto-align** before manual adjustments
2. **Lower thresholds** if you see too much black (masked) area
3. **Use Log₂ ratio** for publication figures (symmetric scale)
4. **Export PNG at 200 DPI** for high-quality figures
5. **Check coverage %** to ensure metabolite is detected

---

For detailed documentation, see [README.md](README.md).
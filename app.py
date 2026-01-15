"""
SPARTA: Spatial Metabolite Alignment & Co-localization Analysis
Main Streamlit Application
"""

import streamlit as st
from pyimzml.ImzMLParser import ImzMLParser
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
import io
import os

# Import backend modules
from backend.centroid_engine import calculate_centroid
from backend.shifter import rigid_shift_image
from backend.rgb_composite import create_rgb_composite, create_rgb_composite_with_green_blend
from backend.colocalization import compute_colocalization_map
from backend.heatmap import reconstruct_heatmap

# Page configuration
st.set_page_config(layout="wide", page_title="SPARTA Metabolite Co-localization")

# Title
st.title("🗺️ SPARTA: Spatial Metabolite Alignment & Co-localization Analysis")


def truncate_name(name, max_len=25):
    """Truncate long metabolite names for cleaner display."""
    if name is None:
        return "N/A"
    name = str(name)
    if len(name) > max_len:
        return name[:max_len-3] + "..."
    return name


def load_imzml_file(path, label=""):
    """Helper to load imzML file with error handling."""
    if not os.path.exists(path):
        return None, f"File not found: {path}"
    try:
        parser = ImzMLParser(path)
        if parser is None or len(parser.coordinates) == 0:
            return None, "imzML loaded but contains no coordinates"
        return parser, None
    except IndexError:
        ibd_path = path.replace('.imzML', '.ibd')
        return None, f"Index error - check .ibd file exists at: {ibd_path}"
    except Exception as e:
        return None, str(e)


# Initialize session state
defaults = {
    'shift_x': 0,
    'shift_y': 0,
    'parser': None,
    'parser_2': None,
    'annotations': pd.DataFrame(),
    'annotations_2': pd.DataFrame(),
    'heatmap_cache': {},
    'heatmap_cache_2': {},
    'imzml_path': None,
    'imzml_path_2': None,
    'analysis_mode': 'Single MSI'
}

for key, default in defaults.items():
    if key not in st.session_state:
        st.session_state[key] = default

# ===== SIDEBAR: MODE SELECTION =====
with st.sidebar:
    st.header("⚙️ Analysis Mode")
    analysis_mode = st.radio(
        "Select Mode",
        options=['Single MSI', 'Comparison Mode'],
        index=0 if st.session_state.analysis_mode == 'Single MSI' else 1,
        help="Single MSI: Analyze metabolite co-localization within one image.\n"
             "Comparison Mode: Compare same metabolite across two MSIs (e.g., pre/post tumor)."
    )
    st.session_state.analysis_mode = analysis_mode
    
    if analysis_mode == 'Comparison Mode':
        st.info("📊 **Comparison Mode**: Load two MSIs to compare the same metabolite between conditions (e.g., pre-tumor vs post-tumor)")
    
    st.divider()

# ===== SIDEBAR: DATA LOADING =====
with st.sidebar:
    if analysis_mode == 'Single MSI':
        st.header("📥 Data Upload")
    else:
        st.header("📥 MSI 1 (Pre-condition)")
    
    # File uploader for imzML 1
    use_local_path = st.checkbox("Use local file path", key="use_local_1")
    if use_local_path:
        local_imzml_path = st.text_input("Local .imzML path", value="data/mouse1.imzML", key="path_1")
        if st.button("Load MSI 1", key="load_1"):
            parser, error = load_imzml_file(local_imzml_path)
            if error:
                st.error(f"Error: {error}")
            else:
                st.session_state.parser = parser
                st.session_state.imzml_path = local_imzml_path
                st.session_state.heatmap_cache = {}
                st.success(f"✅ Loaded {len(parser.coordinates)} pixels")
    else:
        imzml_file = st.file_uploader("Upload .imzML file", type=["imzml", "ibd"], key="upload_1")
        if imzml_file is not None:
            try:
                with open("temp_upload_1.imzML", "wb") as f:
                    f.write(imzml_file.getbuffer())
                parser, error = load_imzml_file("temp_upload_1.imzML")
                if error:
                    st.error(f"Error: {error}")
                else:
                    st.session_state.parser = parser
                    st.session_state.imzml_path = "temp_upload_1.imzML"
                    st.session_state.heatmap_cache = {}
                    st.success(f"✅ Loaded {len(parser.coordinates)} pixels")
            except Exception as e:
                st.error(f"Error: {e}")
    
    # CSV 1
    csv_label = "Upload METASPACE CSV" if analysis_mode == 'Single MSI' else "Upload CSV (MSI 1)"
    csv_file = st.file_uploader(csv_label, type="csv", key="csv_1")
    if csv_file is not None:
        try:
            annotations = pd.read_csv(csv_file, sep=",", engine="python", comment="#")
            annotations = annotations.rename(columns={"moleculeNames": "name"})
            st.session_state.annotations = annotations
            st.success(f"✅ Loaded {len(annotations)} metabolites")
        except Exception as e:
            st.error(f"Error loading CSV: {e}")

# ===== SIDEBAR: SECOND MSI LOADING (Comparison Mode) =====
if analysis_mode == 'Comparison Mode':
    with st.sidebar:
        st.divider()
        st.header("📥 MSI 2 (Post-condition)")
        
        use_local_path_2 = st.checkbox("Use local file path", key="use_local_2")
        if use_local_path_2:
            local_imzml_path_2 = st.text_input("Local .imzML path", value="data/mouse2.imzML", key="path_2")
            if st.button("Load MSI 2", key="load_2"):
                parser_2, error = load_imzml_file(local_imzml_path_2)
                if error:
                    st.error(f"Error: {error}")
                else:
                    st.session_state.parser_2 = parser_2
                    st.session_state.imzml_path_2 = local_imzml_path_2
                    st.session_state.heatmap_cache_2 = {}
                    st.success(f"✅ Loaded {len(parser_2.coordinates)} pixels")
        else:
            imzml_file_2 = st.file_uploader("Upload .imzML file", type=["imzml", "ibd"], key="upload_2")
            if imzml_file_2 is not None:
                try:
                    with open("temp_upload_2.imzML", "wb") as f:
                        f.write(imzml_file_2.getbuffer())
                    parser_2, error = load_imzml_file("temp_upload_2.imzML")
                    if error:
                        st.error(f"Error: {error}")
                    else:
                        st.session_state.parser_2 = parser_2
                        st.session_state.imzml_path_2 = "temp_upload_2.imzML"
                        st.session_state.heatmap_cache_2 = {}
                        st.success(f"✅ Loaded {len(parser_2.coordinates)} pixels")
                except Exception as e:
                    st.error(f"Error: {e}")
        
        # CSV 2
        csv_file_2 = st.file_uploader("Upload CSV (MSI 2)", type="csv", key="csv_2")
        if csv_file_2 is not None:
            try:
                annotations_2 = pd.read_csv(csv_file_2, sep=",", engine="python", comment="#")
                annotations_2 = annotations_2.rename(columns={"moleculeNames": "name"})
                st.session_state.annotations_2 = annotations_2
                st.success(f"✅ Loaded {len(annotations_2)} metabolites")
            except Exception as e:
                st.error(f"Error loading CSV: {e}")

# Check if data is loaded
if st.session_state.parser is None:
    st.warning("⚠️ Please upload an .imzML file to begin.")
    st.info("💡 Tip: Use the local file path option to load `data/mouse1.imzML`")
    st.stop()

if analysis_mode == 'Comparison Mode' and st.session_state.parser_2 is None:
    st.warning("⚠️ Comparison Mode requires two MSI files. Please load MSI 2.")
    st.stop()

parser = st.session_state.parser
parser_2 = st.session_state.parser_2
annotations = st.session_state.annotations
annotations_2 = st.session_state.annotations_2

# ===============================
# SINGLE MSI MODE
# ===============================
if analysis_mode == 'Single MSI':
    st.markdown("**Mode**: Single MSI Co-localization Analysis")
    
    # ===== SIDEBAR: METABOLITE SELECTION =====
    with st.sidebar:
        st.divider()
        st.header("🔬 Metabolite Selection")
        
        # Build metabolite list
        if not annotations.empty and 'name' in annotations.columns:
            metabolite_list = annotations['name'].tolist()
            has_annotations = True
        else:
            metabolite_list = None
            has_annotations = False
        
        # Metabolite A
        st.subheader("Metabolite A")
        if has_annotations and metabolite_list:
            search_a = st.selectbox("Select Metabolite A", metabolite_list, key="metabolite_a")
            row_a = annotations[annotations['name'] == search_a].iloc[0].to_dict()
        else:
            mz_a = st.number_input("m/z A", value=500.0, step=0.1, key="mz_a")
            row_a = {'mz': mz_a, 'name': f"m/z {mz_a:.2f}", 'formula': 'Unknown'}
        
        if 'mz' in row_a:
            st.metric("m/z", f"{row_a['mz']:.4f}")
        st.caption(f"**{truncate_name(row_a.get('name', 'N/A'), 40)}**")
        
        # Compute heatmap A
        mz_a_key = row_a.get('mz', 500.0)
        if mz_a_key not in st.session_state.heatmap_cache:
            with st.spinner("Computing..."):
                try:
                    heatmap_a = reconstruct_heatmap(parser, mz_a_key, tolerance=0.5)
                    st.session_state.heatmap_cache[mz_a_key] = heatmap_a
                except Exception as e:
                    st.error(f"Error: {e}")
                    st.stop()
        else:
            heatmap_a = st.session_state.heatmap_cache[mz_a_key]
        
        coverage_a = (heatmap_a > 0).sum() / heatmap_a.size * 100
        mean_int_a = heatmap_a[heatmap_a > 0].mean() if (heatmap_a > 0).any() else 0
        st.metric("Coverage", f"{coverage_a:.1f}%")
        
        st.divider()
        
        # Metabolite B
        st.subheader("Metabolite B")
        if has_annotations and metabolite_list and len(metabolite_list) > 1:
            default_idx = min(1, len(metabolite_list) - 1)
            search_b = st.selectbox("Select Metabolite B", metabolite_list, key="metabolite_b", index=default_idx)
            row_b = annotations[annotations['name'] == search_b].iloc[0].to_dict()
        else:
            mz_b = st.number_input("m/z B", value=600.0, step=0.1, key="mz_b")
            row_b = {'mz': mz_b, 'name': f"m/z {mz_b:.2f}", 'formula': 'Unknown'}
        
        if 'mz' in row_b:
            st.metric("m/z", f"{row_b['mz']:.4f}")
        st.caption(f"**{truncate_name(row_b.get('name', 'N/A'), 40)}**")
        
        # Compute heatmap B
        mz_b_key = row_b.get('mz', 600.0)
        if mz_b_key not in st.session_state.heatmap_cache:
            with st.spinner("Computing..."):
                try:
                    heatmap_b = reconstruct_heatmap(parser, mz_b_key, tolerance=0.5)
                    st.session_state.heatmap_cache[mz_b_key] = heatmap_b
                except Exception as e:
                    st.error(f"Error: {e}")
                    st.stop()
        else:
            heatmap_b = st.session_state.heatmap_cache[mz_b_key]
        
        coverage_b = (heatmap_b > 0).sum() / heatmap_b.size * 100
        st.metric("Coverage", f"{coverage_b:.1f}%")

    # ===== SIDEBAR: ALIGNMENT CONTROLS =====
    with st.sidebar:
        st.divider()
        st.header("🎯 Alignment Controls")
        
        centroid_method = st.radio(
            "Centroid method",
            options=['Intensity-Weighted', 'Binary-Mask'],
            help="Weighted: emphasizes high-intensity regions. Binary: geometric center."
        )
        centroid_method_key = 'weighted' if centroid_method == 'Intensity-Weighted' else 'binary'
        
        snr_multiplier = st.slider("SNR Floor (k)", 1.0, 5.0, 3.0, 0.5)
        
        if st.button("🔄 Auto-Align"):
            with st.spinner("Calculating centroids..."):
                cx_a, cy_a, meta_a = calculate_centroid(heatmap_a, method=centroid_method_key, snr_floor_multiplier=snr_multiplier)
                cx_b, cy_b, meta_b = calculate_centroid(heatmap_b, method=centroid_method_key, snr_floor_multiplier=snr_multiplier)
            
            if cx_a is not None and cx_b is not None:
                st.session_state.shift_x = int(round(cx_b - cx_a))
                st.session_state.shift_y = int(round(cy_b - cy_a))
                st.success(f"✅ Δx={st.session_state.shift_x}, Δy={st.session_state.shift_y}")
            else:
                st.error("❌ Centroid calculation failed")
        
        # Manual nudge
        st.write("**Manual Nudge**")
        c1, c2, c3, c4 = st.columns(4)
        with c1:
            if st.button("⬅️"):
                st.session_state.shift_x -= 1
        with c2:
            if st.button("⬆️"):
                st.session_state.shift_y -= 1
        with c3:
            if st.button("⬇️"):
                st.session_state.shift_y += 1
        with c4:
            if st.button("➡️"):
                st.session_state.shift_x += 1
        
        st.metric("Δx", st.session_state.shift_x)
        st.metric("Δy", st.session_state.shift_y)
        
        if st.button("🔁 Reset"):
            st.session_state.shift_x = 0
            st.session_state.shift_y = 0

    # ===== MAIN: ALIGNMENT VALIDATION =====
    st.header("🔍 Alignment Validation")
    
    try:
        heatmap_b_shifted = rigid_shift_image(heatmap_b, st.session_state.shift_x, st.session_state.shift_y)
    except ValueError as e:
        st.error(f"Shift error: {e}")
        heatmap_b_shifted = heatmap_b.copy()
    
    # Cleaner titles
    title_a = truncate_name(row_a.get('name', 'A'), 20)
    title_b = truncate_name(row_b.get('name', 'B'), 20)
    
    fig_val = make_subplots(
        rows=1, cols=3,
        subplot_titles=[f"A: {title_a}", f"B (Aligned): {title_b}", "Overlay (R=A, B=B)"],
        specs=[[{"type": "heatmap"}, {"type": "heatmap"}, {"type": "image"}]]
    )
    
    fig_val.add_trace(go.Heatmap(z=heatmap_a, colorscale='Viridis', showscale=False), row=1, col=1)
    fig_val.add_trace(go.Heatmap(z=heatmap_b_shifted, colorscale='Viridis', showscale=False), row=1, col=2)
    
    overlay = create_rgb_composite(heatmap_a, heatmap_b_shifted, percentile_clip=99.0)
    overlay_uint8 = (overlay * 255).astype(np.uint8)
    fig_val.add_trace(go.Image(z=overlay_uint8), row=1, col=3)
    
    fig_val.update_layout(height=400, showlegend=False, margin=dict(t=50, b=30))
    st.plotly_chart(fig_val, use_container_width=True)

    # ===== MAIN: CO-LOCALIZATION MAP =====
    st.header("🗺️ Co-localization Map")
    
    col1, col2 = st.columns(2)
    with col1:
        min_b_pct = st.slider("Min B threshold (%)", 0.0, 100.0, 10.0, 1.0)
        max_b_int = float(heatmap_b_shifted.max())
        min_b_int = (min_b_pct / 100.0) * max_b_int if max_b_int > 0 else 0.0
    with col2:
        handle_zeros = st.radio("Masked pixels", options=['Transparent', 'White'])
    
    colocalization_map = compute_colocalization_map(
        heatmap_a, heatmap_b_shifted,
        min_b_intensity=min_b_int,
        handle_zeros='zero' if handle_zeros == 'White' else 'nan'
    )
    
    valid_mask = ~np.isnan(colocalization_map)
    if np.any(valid_mask):
        coloc_log = np.full_like(colocalization_map, np.nan)
        coloc_log[valid_mask] = np.log2(np.clip(colocalization_map[valid_mask], 0.01, 100))
        
        fig_coloc = go.Figure(data=go.Heatmap(
            z=coloc_log,
            colorscale='RdBu_r',
            zmid=0,
            colorbar=dict(title="log₂(A/B)", tickvals=[-2, -1, 0, 1, 2], ticktext=['0.25', '0.5', '1', '2', '4'])
        ))
        
        fig_coloc.update_layout(
            title=f"{title_a} vs {title_b}",
            xaxis_title="X (px)", yaxis_title="Y (px)",
            height=550,
            yaxis=dict(scaleanchor="x"),
            margin=dict(t=50, b=30)
        )
        st.plotly_chart(fig_coloc, use_container_width=True)
        
        st.info("🔴 Red = A dominant | 🔵 Blue = B dominant | ⚪ White = Balanced")
    else:
        st.warning("❌ No valid co-localization regions.")

    # ===== SIDEBAR: EXPORT =====
    with st.sidebar:
        st.divider()
        st.header("💾 Export")
        
        if st.button("📊 Export CSV"):
            y_coords, x_coords = np.indices(colocalization_map.shape)
            export_df = pd.DataFrame({
                'x': x_coords.flatten(),
                'y': y_coords.flatten(),
                f'{row_a.get("name", "A")[:20]}_int': heatmap_a.flatten(),
                f'{row_b.get("name", "B")[:20]}_int': heatmap_b_shifted.flatten(),
                'coloc_index': colocalization_map.flatten()
            })
            csv = export_df.to_csv(index=False)
            st.download_button("Download", csv, "SPARTA_colocalization.csv", "text/csv")
        
        if st.button("🖼️ Export PNG"):
            fig, axes = plt.subplots(1, 4, figsize=(16, 4))
            
            axes[0].imshow(heatmap_a, cmap='viridis')
            axes[0].set_title(f"A: {title_a}", fontsize=10)
            axes[0].axis('off')
            
            axes[1].imshow(heatmap_b_shifted, cmap='viridis')
            axes[1].set_title(f"B: {title_b}", fontsize=10)
            axes[1].axis('off')
            
            axes[2].imshow(overlay)
            axes[2].set_title("Overlay", fontsize=10)
            axes[2].axis('off')
            
            if np.any(valid_mask):
                im = axes[3].imshow(coloc_log, cmap='RdBu_r', vmin=-2, vmax=2)
                axes[3].set_title("Co-loc (log₂)", fontsize=10)
                plt.colorbar(im, ax=axes[3], fraction=0.046)
            else:
                axes[3].text(0.5, 0.5, 'No data', ha='center', va='center')
            axes[3].axis('off')
            
            plt.tight_layout()
            buf = io.BytesIO()
            plt.savefig(buf, format='png', dpi=200, bbox_inches='tight')
            buf.seek(0)
            st.download_button("Download", buf.getvalue(), "SPARTA_analysis.png", "image/png")
            st.success("✅ PNG ready!")

# ===============================
# COMPARISON MODE
# ===============================
else:
    st.markdown("**Mode**: MSI Comparison (Pre vs Post Condition)")
    
    # ===== SIDEBAR: METABOLITE SELECTION (COMPARISON) =====
    with st.sidebar:
        st.divider()
        st.header("🔬 Metabolite Selection")
        
        # Build metabolite lists from both annotations
        has_ann_1 = not annotations.empty and 'name' in annotations.columns
        has_ann_2 = not annotations_2.empty and 'name' in annotations_2.columns
        
        if has_ann_1:
            mlist_1 = annotations['name'].tolist()
        else:
            mlist_1 = None
        
        if has_ann_2:
            mlist_2 = annotations_2['name'].tolist()
        else:
            mlist_2 = None
        
        # Find common metabolites if both have annotations
        if mlist_1 and mlist_2:
            common_metabolites = sorted(set(mlist_1) & set(mlist_2))
            st.info(f"📋 {len(common_metabolites)} common metabolites found")
            use_common = st.checkbox("Show only common metabolites", value=True)
            
            if use_common and common_metabolites:
                target_mlist = common_metabolites
            else:
                target_mlist = mlist_1
        else:
            target_mlist = mlist_1
            use_common = False
        
        st.subheader("Target Metabolite")
        
        if target_mlist:
            selected_metabolite = st.selectbox("Select Metabolite", target_mlist, key="comp_metabolite")
            
            # Get row from annotations 1
            if has_ann_1 and selected_metabolite in annotations['name'].values:
                row_1 = annotations[annotations['name'] == selected_metabolite].iloc[0].to_dict()
            else:
                row_1 = {'name': selected_metabolite, 'mz': 500.0}
            
            # Get row from annotations 2 (may have different m/z)
            if has_ann_2 and selected_metabolite in annotations_2['name'].values:
                row_2 = annotations_2[annotations_2['name'] == selected_metabolite].iloc[0].to_dict()
            else:
                row_2 = row_1.copy()
        else:
            target_mz = st.number_input("Target m/z", value=500.0, step=0.1, key="comp_mz")
            row_1 = {'mz': target_mz, 'name': f"m/z {target_mz:.2f}"}
            row_2 = row_1.copy()
            selected_metabolite = row_1['name']
        
        mz_1 = row_1.get('mz', 500.0)
        mz_2 = row_2.get('mz', mz_1)
        
        st.metric("m/z (MSI 1)", f"{mz_1:.4f}")
        st.metric("m/z (MSI 2)", f"{mz_2:.4f}")
        st.caption(f"**{truncate_name(selected_metabolite, 35)}**")

    # Compute heatmaps for both MSIs
    cache_key_1 = (1, mz_1)
    cache_key_2 = (2, mz_2)
    
    if cache_key_1 not in st.session_state.heatmap_cache:
        with st.spinner("Computing heatmap for MSI 1..."):
            heatmap_1 = reconstruct_heatmap(parser, mz_1, tolerance=0.5)
            st.session_state.heatmap_cache[cache_key_1] = heatmap_1
    else:
        heatmap_1 = st.session_state.heatmap_cache[cache_key_1]
    
    if cache_key_2 not in st.session_state.heatmap_cache_2:
        with st.spinner("Computing heatmap for MSI 2..."):
            heatmap_2 = reconstruct_heatmap(parser_2, mz_2, tolerance=0.5)
            st.session_state.heatmap_cache_2[cache_key_2] = heatmap_2
    else:
        heatmap_2 = st.session_state.heatmap_cache_2[cache_key_2]

    # ===== SIDEBAR: ALIGNMENT CONTROLS (COMPARISON) =====
    with st.sidebar:
        st.divider()
        st.header("🎯 Alignment Controls")
        
        centroid_method = st.radio(
            "Centroid method",
            options=['Intensity-Weighted', 'Binary-Mask'],
            key="comp_centroid"
        )
        centroid_method_key = 'weighted' if centroid_method == 'Intensity-Weighted' else 'binary'
        
        snr_mult = st.slider("SNR Floor (k)", 1.0, 5.0, 3.0, 0.5, key="comp_snr")
        
        if st.button("🔄 Auto-Align MSIs"):
            with st.spinner("Calculating centroids..."):
                cx_1, cy_1, meta_1 = calculate_centroid(heatmap_1, method=centroid_method_key, snr_floor_multiplier=snr_mult)
                cx_2, cy_2, meta_2 = calculate_centroid(heatmap_2, method=centroid_method_key, snr_floor_multiplier=snr_mult)
            
            if cx_1 is not None and cx_2 is not None:
                st.session_state.shift_x = int(round(cx_2 - cx_1))
                st.session_state.shift_y = int(round(cy_2 - cy_1))
                st.success(f"✅ Δx={st.session_state.shift_x}, Δy={st.session_state.shift_y}")
            else:
                st.error("❌ Centroid calculation failed")
        
        st.write("**Manual Nudge**")
        c1, c2, c3, c4 = st.columns(4)
        with c1:
            if st.button("⬅️", key="comp_l"):
                st.session_state.shift_x -= 1
        with c2:
            if st.button("⬆️", key="comp_u"):
                st.session_state.shift_y -= 1
        with c3:
            if st.button("⬇️", key="comp_d"):
                st.session_state.shift_y += 1
        with c4:
            if st.button("➡️", key="comp_r"):
                st.session_state.shift_x += 1
        
        st.metric("Δx", st.session_state.shift_x)
        st.metric("Δy", st.session_state.shift_y)
        
        if st.button("🔁 Reset", key="comp_reset"):
            st.session_state.shift_x = 0
            st.session_state.shift_y = 0

    # Resize heatmaps to match if needed
    h1, w1 = heatmap_1.shape
    h2, w2 = heatmap_2.shape
    max_h, max_w = max(h1, h2), max(w1, w2)
    
    # Pad smaller heatmap to match larger
    hm1_padded = np.zeros((max_h, max_w), dtype=np.float32)
    hm2_padded = np.zeros((max_h, max_w), dtype=np.float32)
    hm1_padded[:h1, :w1] = heatmap_1
    hm2_padded[:h2, :w2] = heatmap_2
    
    # Apply shift to MSI 2
    try:
        hm2_shifted = rigid_shift_image(hm2_padded, st.session_state.shift_x, st.session_state.shift_y)
    except ValueError:
        hm2_shifted = hm2_padded.copy()

    # ===== MAIN: COMPARISON PANEL =====
    st.header("📊 Comparison Analysis Panel")
    
    metabolite_short = truncate_name(selected_metabolite, 25)
    
    # Row 1: Side-by-side heatmaps
    st.subheader("Individual Heatmaps")
    
    fig_side = make_subplots(
        rows=1, cols=2,
        subplot_titles=["MSI 1 (Pre)", "MSI 2 (Post, Aligned)"],
        horizontal_spacing=0.08
    )
    
    fig_side.add_trace(go.Heatmap(z=hm1_padded, colorscale='Viridis', showscale=True, name="Pre"), row=1, col=1)
    fig_side.add_trace(go.Heatmap(z=hm2_shifted, colorscale='Viridis', showscale=True, name="Post"), row=1, col=2)
    
    fig_side.update_layout(
        height=400,
        title=f"{metabolite_short}",
        margin=dict(t=60, b=30)
    )
    st.plotly_chart(fig_side, use_container_width=True)

    # Row 2: Ratio Heatmap & Overlay
    st.subheader("Ratio Analysis")
    
    col1, col2 = st.columns(2)
    with col1:
        ratio_type = st.radio("Ratio Type", options=['Post/Pre', 'Pre/Post', 'Log₂(Post/Pre)'], horizontal=True)
    with col2:
        min_thresh_pct = st.slider("Min intensity threshold (%)", 0.0, 50.0, 5.0, 1.0, key="comp_thresh")
    
    # Compute threshold
    combined_max = max(hm1_padded.max(), hm2_shifted.max())
    min_thresh = (min_thresh_pct / 100.0) * combined_max
    
    # Create valid mask (both must have signal)
    valid_mask = (hm1_padded > min_thresh) & (hm2_shifted > min_thresh)
    
    # Compute ratio map
    ratio_map = np.full_like(hm1_padded, np.nan)
    
    if ratio_type == 'Post/Pre':
        ratio_map[valid_mask] = hm2_shifted[valid_mask] / (hm1_padded[valid_mask] + 1e-9)
        cbar_title = "Post/Pre"
        colorscale = 'RdBu_r'
        zmid = 1.0
    elif ratio_type == 'Pre/Post':
        ratio_map[valid_mask] = hm1_padded[valid_mask] / (hm2_shifted[valid_mask] + 1e-9)
        cbar_title = "Pre/Post"
        colorscale = 'RdBu_r'
        zmid = 1.0
    else:  # Log2
        ratio_raw = hm2_shifted[valid_mask] / (hm1_padded[valid_mask] + 1e-9)
        ratio_map[valid_mask] = np.log2(np.clip(ratio_raw, 0.01, 100))
        cbar_title = "log₂(Post/Pre)"
        colorscale = 'RdBu_r'
        zmid = 0.0
    
    # Create panel
    fig_ratio = make_subplots(
        rows=1, cols=3,
        subplot_titles=["Ratio Heatmap", "Overlay (R=Pre, B=Post)", "Difference Map"],
        specs=[[{"type": "heatmap"}, {"type": "image"}, {"type": "heatmap"}]],
        horizontal_spacing=0.06
    )
    
    # Ratio heatmap
    if ratio_type == 'Log₂(Post/Pre)':
        fig_ratio.add_trace(go.Heatmap(
            z=ratio_map, colorscale=colorscale, zmid=zmid, showscale=True,
            colorbar=dict(title=cbar_title, x=0.30)
        ), row=1, col=1)
    else:
        # Clip for visualization
        ratio_clipped = np.clip(ratio_map, 0.1, 10)
        fig_ratio.add_trace(go.Heatmap(
            z=ratio_clipped, colorscale=colorscale, zmid=zmid, showscale=True,
            colorbar=dict(title=cbar_title, x=0.30)
        ), row=1, col=1)
    
    # RGB overlay
    overlay_comp = create_rgb_composite(hm1_padded, hm2_shifted, percentile_clip=99.0)
    overlay_uint8 = (overlay_comp * 255).astype(np.uint8)
    fig_ratio.add_trace(go.Image(z=overlay_uint8), row=1, col=2)
    
    # Difference map (Post - Pre, normalized)
    diff_map = np.full_like(hm1_padded, np.nan)
    if combined_max > 0:
        diff_map[valid_mask] = (hm2_shifted[valid_mask] - hm1_padded[valid_mask]) / combined_max
    fig_ratio.add_trace(go.Heatmap(
        z=diff_map, colorscale='RdBu_r', zmid=0, showscale=True,
        colorbar=dict(title="Δ (norm)", x=1.0)
    ), row=1, col=3)
    
    fig_ratio.update_layout(
        height=450,
        margin=dict(t=60, b=30)
    )
    st.plotly_chart(fig_ratio, use_container_width=True)
    
    st.info("""
    **🔴 Red regions**: Higher in Pre (MSI 1) | **🔵 Blue regions**: Higher in Post (MSI 2) | **⚪ White**: Similar intensity
    """)

    # Statistics
    st.subheader("📈 Comparison Statistics")
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Pre Mean", f"{hm1_padded[hm1_padded > 0].mean():.2e}" if (hm1_padded > 0).any() else "0")
        st.metric("Pre Max", f"{hm1_padded.max():.2e}")
    
    with col2:
        st.metric("Post Mean", f"{hm2_shifted[hm2_shifted > 0].mean():.2e}" if (hm2_shifted > 0).any() else "0")
        st.metric("Post Max", f"{hm2_shifted.max():.2e}")
    
    with col3:
        cov_1 = (hm1_padded > min_thresh).sum() / hm1_padded.size * 100
        cov_2 = (hm2_shifted > min_thresh).sum() / hm2_shifted.size * 100
        st.metric("Pre Coverage", f"{cov_1:.1f}%")
        st.metric("Post Coverage", f"{cov_2:.1f}%")
    
    with col4:
        if np.any(valid_mask):
            if ratio_type == 'Log₂(Post/Pre)':
                mean_ratio = np.nanmean(ratio_map[valid_mask])
                st.metric("Mean log₂ Ratio", f"{mean_ratio:.2f}")
                fold_change = 2 ** mean_ratio
                st.metric("Fold Change", f"{fold_change:.2f}x")
            else:
                mean_ratio = np.nanmean(ratio_map[valid_mask])
                st.metric("Mean Ratio", f"{mean_ratio:.2f}")
                st.metric("Valid Pixels", f"{valid_mask.sum()}")
        else:
            st.metric("Mean Ratio", "N/A")
            st.metric("Valid Pixels", "0")

    # Centroid visualization
    st.subheader("🎯 Centroid Mapping")
    
    cx_1, cy_1, meta_1 = calculate_centroid(hm1_padded, method=centroid_method_key, snr_floor_multiplier=snr_mult)
    cx_2, cy_2, meta_2 = calculate_centroid(hm2_shifted, method=centroid_method_key, snr_floor_multiplier=snr_mult)
    
    fig_centroid = make_subplots(
        rows=1, cols=2,
        subplot_titles=["Pre + Centroid", "Post + Centroid"],
        horizontal_spacing=0.08
    )
    
    fig_centroid.add_trace(go.Heatmap(z=hm1_padded, colorscale='Viridis', showscale=False), row=1, col=1)
    fig_centroid.add_trace(go.Heatmap(z=hm2_shifted, colorscale='Viridis', showscale=False), row=1, col=2)
    
    # Add centroid markers
    if cx_1 is not None:
        fig_centroid.add_trace(go.Scatter(
            x=[cx_1], y=[cy_1], mode='markers',
            marker=dict(size=15, color='red', symbol='x', line=dict(width=2, color='white')),
            name='Pre Centroid'
        ), row=1, col=1)
    
    if cx_2 is not None:
        fig_centroid.add_trace(go.Scatter(
            x=[cx_2], y=[cy_2], mode='markers',
            marker=dict(size=15, color='cyan', symbol='x', line=dict(width=2, color='white')),
            name='Post Centroid'
        ), row=1, col=2)
    
    fig_centroid.update_layout(height=400, showlegend=True, margin=dict(t=60, b=30))
    st.plotly_chart(fig_centroid, use_container_width=True)
    
    if cx_1 is not None and cx_2 is not None:
        dist = np.sqrt((cx_2 - cx_1)**2 + (cy_2 - cy_1)**2)
        st.info(f"📍 **Centroid shift**: Δx={cx_2-cx_1:.1f}, Δy={cy_2-cy_1:.1f}, Distance={dist:.1f} px")

    # ===== SIDEBAR: EXPORT (COMPARISON) =====
    with st.sidebar:
        st.divider()
        st.header("💾 Export")
        
        if st.button("📊 Export Comparison CSV"):
            y_coords, x_coords = np.indices(ratio_map.shape)
            export_df = pd.DataFrame({
                'x': x_coords.flatten(),
                'y': y_coords.flatten(),
                'pre_intensity': hm1_padded.flatten(),
                'post_intensity': hm2_shifted.flatten(),
                'ratio': ratio_map.flatten()
            })
            csv = export_df.to_csv(index=False)
            st.download_button("Download", csv, f"SPARTA_comparison_{metabolite_short[:15]}.csv", "text/csv")
        
        if st.button("🖼️ Export Panel PNG"):
            fig, axes = plt.subplots(2, 3, figsize=(14, 9))
            
            axes[0, 0].imshow(hm1_padded, cmap='viridis')
            axes[0, 0].set_title("Pre (MSI 1)", fontsize=10)
            if cx_1:
                axes[0, 0].plot(cx_1, cy_1, 'rx', markersize=10, markeredgewidth=2)
            axes[0, 0].axis('off')
            
            axes[0, 1].imshow(hm2_shifted, cmap='viridis')
            axes[0, 1].set_title("Post (MSI 2)", fontsize=10)
            if cx_2:
                axes[0, 1].plot(cx_2, cy_2, 'cx', markersize=10, markeredgewidth=2)
            axes[0, 1].axis('off')
            
            axes[0, 2].imshow(overlay_comp)
            axes[0, 2].set_title("Overlay", fontsize=10)
            axes[0, 2].axis('off')
            
            if ratio_type == 'Log₂(Post/Pre)':
                im = axes[1, 0].imshow(ratio_map, cmap='RdBu_r', vmin=-3, vmax=3)
            else:
                im = axes[1, 0].imshow(np.clip(ratio_map, 0.1, 10), cmap='RdBu_r', vmin=0.1, vmax=10)
            axes[1, 0].set_title(f"Ratio ({cbar_title})", fontsize=10)
            plt.colorbar(im, ax=axes[1, 0], fraction=0.046)
            axes[1, 0].axis('off')
            
            im2 = axes[1, 1].imshow(diff_map, cmap='RdBu_r', vmin=-0.5, vmax=0.5)
            axes[1, 1].set_title("Difference (norm)", fontsize=10)
            plt.colorbar(im2, ax=axes[1, 1], fraction=0.046)
            axes[1, 1].axis('off')
            
            # Histogram
            if np.any(valid_mask):
                valid_ratios = ratio_map[valid_mask]
                valid_ratios = valid_ratios[~np.isnan(valid_ratios)]
                if len(valid_ratios) > 0:
                    axes[1, 2].hist(valid_ratios, bins=50, color='steelblue', edgecolor='white')
                    axes[1, 2].axvline(x=np.median(valid_ratios), color='red', linestyle='--', label=f'Median: {np.median(valid_ratios):.2f}')
                    axes[1, 2].set_xlabel(cbar_title)
                    axes[1, 2].set_ylabel("Frequency")
                    axes[1, 2].legend()
            axes[1, 2].set_title("Ratio Distribution", fontsize=10)
            
            plt.suptitle(f"{metabolite_short}", fontsize=12, fontweight='bold')
            plt.tight_layout()
            
            buf = io.BytesIO()
            plt.savefig(buf, format='png', dpi=200, bbox_inches='tight')
            buf.seek(0)
            st.download_button("Download", buf.getvalue(), f"SPARTA_comparison_{metabolite_short[:15]}.png", "image/png")
            st.success("✅ PNG ready!")

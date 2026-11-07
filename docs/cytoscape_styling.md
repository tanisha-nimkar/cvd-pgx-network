# Cytoscape Network Visualization & Styling Protocol

**Source Graph:** `results/network/network.graphml`  
**Network Attributes:** 135 nodes, 225 edges (STRING v12.5 text-mining-excluded functional LCC)  
**Checksum:** `33bcd28abc5583d264ebda841b09cdba2364ceeccd5d4d63b5a6501b068c8b23`

---

## Visual Encoding Scheme

| Visual Property | Target Element | Mapping Attribute | Visual Style / Value |
|:---|:---|:---|:---|
| **Node Shape** | Panel Seeds | `is_panel_gene == True` | Diamond / Octagon |
| **Node Shape** | Interactors | `is_panel_gene == False` | Ellipse |
| **Node Size** | All Nodes | `degree` | Continuous: 25 px (degree=1) to 75 px (degree=89) |
| **Node Fill Color**| Panel Seeds | `gene_symbol` | `#d95f02` (Accent Orange) |
| **Node Fill Color**| 1st-Degree Interactors | Topological Hubs | `#2b5c8f` (Navy Blue) |
| **Node Label** | All Nodes | `gene_symbol` | Font: SansSerif, 10pt (Bold for panel seeds) |
| **Edge Width** | All Edges | `combined_score` | Continuous: 1.0 (score=0.70) to 3.5 (score=0.99) |
| **Edge Color** | All Edges | Default | `#999999` with 60% opacity |
| **Layout Algorithm**| All Nodes | Prefuse Force Directed | Spring length: 50, Spring strength: 0.0001 |

---

## Session File & Export Specifications
- **Cytoscape Version:** >= 3.10.0
- **Export Resolution:** 300 DPI, lossless PNG / PDF vector export
- **Export Target:** `figures/fig4_cytoscape_network.png`
- **Topology Rule:** Absolute freeze — zero edge pruning, zero manual node additions.

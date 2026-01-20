#!/usr/bin/env python3
"""
Quick Three-Solver Comparison: XFOIL vs Mfoil vs RustFoil

Uses existing XFOIL/RustFoil data and adds Mfoil results.
"""

import json
import sys
import warnings
from pathlib import Path
import numpy as np

# Add mfoil to path
PROJECT_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_DIR / "mfoil"))

warnings.filterwarnings('ignore')
from mfoil import mfoil

OUTPUT_DIR = PROJECT_DIR / "testdata" / "xfoil_comparison"

# Test configurations (subset for speed)
CONFIGS = [
    ("0012", 1e6),
    ("0012", 3e6),
    ("2412", 1e6),
    ("2412", 3e6),
    ("4412", 1e6),
    ("4412", 3e6),
]

ALPHAS = list(range(-4, 13, 1))  # Up to 12 degrees to avoid separation issues


def run_mfoil(naca: str, reynolds: float, alphas: list) -> dict:
    """Run Mfoil for given configuration."""
    data = {"alpha": [], "cl": [], "cd": [], "cm": []}
    
    print(f"    Running Mfoil for NACA {naca} Re={reynolds:.0e}...", end=" ", flush=True)
    
    try:
        m = mfoil(naca=naca, npanel=129)
        m.param.verb = 0
        m.param.doplot = False
        m.param.niglob = 20
        
        for alpha in alphas:
            try:
                m.setoper(alpha=alpha, Re=reynolds, Ma=0.0)
                m.oper.initbl = True
                m.solve()
                
                if m.glob.conv and hasattr(m.post, 'cl'):
                    cl = float(m.post.cl)
                    cd = float(m.post.cd) if hasattr(m.post, 'cd') else 0.0
                    cm = float(m.post.cm) if hasattr(m.post, 'cm') else 0.0
                    
                    if np.isfinite(cl) and np.isfinite(cd) and abs(cl) < 5 and cd >= 0 and cd < 1:
                        data["alpha"].append(alpha)
                        data["cl"].append(cl)
                        data["cd"].append(cd)
                        data["cm"].append(cm)
            except:
                pass
        
        print(f"✓ ({len(data['alpha'])} points)")
    except Exception as e:
        print(f"✗ ({e})")
        return None
    
    return data if data["alpha"] else None


def load_existing_data():
    """Load existing XFOIL vs RustFoil comparison data."""
    json_path = OUTPUT_DIR / "xfoil_comparison_data.json"
    if not json_path.exists():
        print("No existing data found")
        return {}
    
    with open(json_path, 'r') as f:
        return json.load(f)


def generate_html(all_data: list, output_path: Path):
    """Generate interactive HTML."""
    
    html = """<!DOCTYPE html>
<html>
<head>
    <title>XFOIL vs Mfoil vs RustFoil Comparison</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        * { box-sizing: border-box; }
        body { 
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; 
            margin: 0; 
            padding: 20px; 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }
        .container { max-width: 1600px; margin: 0 auto; }
        h1 { color: white; text-align: center; margin-bottom: 10px; text-shadow: 2px 2px 4px rgba(0,0,0,0.3); }
        .subtitle { color: rgba(255,255,255,0.9); text-align: center; margin-bottom: 20px; font-size: 14px; }
        .card { background: rgba(255,255,255,0.95); padding: 15px 20px; border-radius: 10px; margin: 15px 0; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }
        .selector { display: flex; align-items: center; gap: 20px; flex-wrap: wrap; }
        .selector label { font-weight: 600; color: #333; }
        select { padding: 10px 15px; font-size: 14px; border: 2px solid #ddd; border-radius: 6px; background: white; cursor: pointer; min-width: 150px; }
        select:hover { border-color: #667eea; }
        .plot-container { display: grid; grid-template-columns: repeat(2, 1fr); gap: 20px; margin-top: 20px; }
        .plot { background: white; padding: 15px; border-radius: 10px; box-shadow: 0 4px 15px rgba(0,0,0,0.15); }
        .legend-info { display: flex; justify-content: center; gap: 30px; }
        .legend-item { display: flex; align-items: center; gap: 8px; }
        .legend-color { width: 30px; height: 4px; border-radius: 2px; }
        .xfoil-color { background: #1f77b4; }
        .mfoil-color { background: #2ca02c; }
        .rustfoil-color { background: #ff7f0e; }
        .stats-table { width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 13px; }
        .stats-table th, .stats-table td { padding: 8px 12px; text-align: left; border-bottom: 1px solid #eee; }
        .stats-table th { background: #f5f5f5; font-weight: 600; }
        @media (max-width: 1200px) { .plot-container { grid-template-columns: 1fr; } }
    </style>
</head>
<body>
    <div class="container">
        <h1>Three-Solver Aerodynamic Comparison</h1>
        <p class="subtitle">XFOIL (Fortran) vs Mfoil (Python) vs RustFoil (Rust)</p>
        
        <div class="card legend-info">
            <div class="legend-item"><div class="legend-color xfoil-color"></div><span><strong>XFOIL</strong> - Reference</span></div>
            <div class="legend-item"><div class="legend-color mfoil-color"></div><span><strong>Mfoil</strong> - Python</span></div>
            <div class="legend-item"><div class="legend-color rustfoil-color"></div><span><strong>RustFoil</strong> - This Project</span></div>
        </div>
        
        <div class="card selector">
            <label>Airfoil:</label>
            <select id="airfoilSelect" onchange="updatePlots()">
                <option value="0012">NACA 0012</option>
                <option value="2412">NACA 2412</option>
                <option value="4412">NACA 4412</option>
            </select>
            <label>Reynolds:</label>
            <select id="reSelect" onchange="updatePlots()">
                <option value="1000000">1 Million</option>
                <option value="3000000">3 Million</option>
            </select>
        </div>
        
        <div class="card" id="statsInfo">
            <h3 style="margin:0 0 10px 0">Comparison Statistics</h3>
            <table class="stats-table" id="statsTable"></table>
        </div>
        
        <div class="plot-container">
            <div class="plot" id="clPlot" style="height: 450px;"></div>
            <div class="plot" id="cdPlot" style="height: 450px;"></div>
            <div class="plot" id="polarPlot" style="height: 450px;"></div>
            <div class="plot" id="ldPlot" style="height: 450px;"></div>
        </div>
    </div>
    
    <script>
        const allData = DATA_PLACEHOLDER;
        
        function computeStats(ref, test) {
            if (!ref || !test || !ref.alpha || !test.alpha || ref.alpha.length === 0) return null;
            let clErrs = [], cdErrs = [];
            for (let i = 0; i < ref.alpha.length; i++) {
                const alpha = ref.alpha[i];
                const j = test.alpha.findIndex(a => Math.abs(a - alpha) < 0.5);
                if (j >= 0 && test.cl[j] != null && test.cd[j] != null) {
                    clErrs.push(test.cl[j] - ref.cl[i]);
                    cdErrs.push(test.cd[j] - ref.cd[i]);
                }
            }
            if (clErrs.length === 0) return null;
            const rms = arr => Math.sqrt(arr.reduce((s, x) => s + x*x, 0) / arr.length);
            const maxAbs = arr => Math.max(...arr.map(Math.abs));
            return { clRms: rms(clErrs), cdRms: rms(cdErrs), clMax: maxAbs(clErrs), cdMax: maxAbs(cdErrs), n: clErrs.length };
        }
        
        function updateStats(xfoil, mfoil, rustfoil) {
            const ms = computeStats(xfoil, mfoil);
            const rs = computeStats(xfoil, rustfoil);
            let html = '<tr><th>Metric</th><th>Mfoil vs XFOIL</th><th>RustFoil vs XFOIL</th></tr>';
            if (ms && rs) {
                html += `<tr><td>Cl RMS Error</td><td>${ms.clRms.toFixed(4)}</td><td>${rs.clRms.toFixed(4)}</td></tr>`;
                html += `<tr><td>Cd RMS (counts)</td><td>${(ms.cdRms*10000).toFixed(1)}</td><td>${(rs.cdRms*10000).toFixed(1)}</td></tr>`;
                html += `<tr><td>Points</td><td>${ms.n}</td><td>${rs.n}</td></tr>`;
            } else {
                html += '<tr><td colspan="3">Insufficient data</td></tr>';
            }
            document.getElementById('statsTable').innerHTML = html;
        }
        
        function updatePlots() {
            const airfoil = document.getElementById('airfoilSelect').value;
            const re = parseFloat(document.getElementById('reSelect').value);
            const match = allData.find(d => d.airfoil === airfoil && d.reynolds === re);
            if (!match) return;
            
            const xf = match.xfoil || {alpha:[], cl:[], cd:[], cm:[]};
            const mf = match.mfoil || {alpha:[], cl:[], cd:[], cm:[]};
            const rf = match.rustfoil || {alpha:[], cl:[], cd:[], cm:[]};
            
            updateStats(xf, mf, rf);
            
            const layout = { font: {family: '-apple-system, sans-serif'}, legend: {x:0.02, y:0.98, bgcolor:'rgba(255,255,255,0.8)'}, margin: {t:50,r:30,b:50,l:60} };
            
            Plotly.newPlot('clPlot', [
                {x:xf.alpha, y:xf.cl, mode:'lines+markers', name:'XFOIL', line:{color:'#1f77b4', width:2}, marker:{size:6}},
                {x:mf.alpha, y:mf.cl, mode:'lines+markers', name:'Mfoil', line:{color:'#2ca02c', width:2, dash:'dash'}, marker:{size:6, symbol:'diamond'}},
                {x:rf.alpha, y:rf.cl, mode:'lines+markers', name:'RustFoil', line:{color:'#ff7f0e', width:2, dash:'dot'}, marker:{size:6, symbol:'square'}}
            ], {...layout, title:'<b>Cl vs Alpha</b>', xaxis:{title:'α (deg)', gridcolor:'#eee'}, yaxis:{title:'Cl', gridcolor:'#eee'}});
            
            Plotly.newPlot('cdPlot', [
                {x:xf.alpha, y:xf.cd.map(c=>c*10000), mode:'lines+markers', name:'XFOIL', line:{color:'#1f77b4', width:2}, marker:{size:6}},
                {x:mf.alpha, y:mf.cd.map(c=>c*10000), mode:'lines+markers', name:'Mfoil', line:{color:'#2ca02c', width:2, dash:'dash'}, marker:{size:6, symbol:'diamond'}},
                {x:rf.alpha, y:rf.cd.map(c=>c*10000), mode:'lines+markers', name:'RustFoil', line:{color:'#ff7f0e', width:2, dash:'dot'}, marker:{size:6, symbol:'square'}}
            ], {...layout, title:'<b>Cd vs Alpha</b>', xaxis:{title:'α (deg)', gridcolor:'#eee'}, yaxis:{title:'Cd (counts)', gridcolor:'#eee'}});
            
            Plotly.newPlot('polarPlot', [
                {x:xf.cd.map(c=>c*10000), y:xf.cl, mode:'lines+markers', name:'XFOIL', line:{color:'#1f77b4', width:2}, marker:{size:6}},
                {x:mf.cd.map(c=>c*10000), y:mf.cl, mode:'lines+markers', name:'Mfoil', line:{color:'#2ca02c', width:2, dash:'dash'}, marker:{size:6, symbol:'diamond'}},
                {x:rf.cd.map(c=>c*10000), y:rf.cl, mode:'lines+markers', name:'RustFoil', line:{color:'#ff7f0e', width:2, dash:'dot'}, marker:{size:6, symbol:'square'}}
            ], {...layout, title:'<b>Drag Polar</b>', xaxis:{title:'Cd (counts)', gridcolor:'#eee'}, yaxis:{title:'Cl', gridcolor:'#eee'}, legend:{x:0.98, y:0.02, xanchor:'right'}});
            
            const ld = (cl,cd) => cl.map((c,i) => cd[i]>1e-6 ? c/cd[i] : null);
            Plotly.newPlot('ldPlot', [
                {x:xf.alpha, y:ld(xf.cl,xf.cd), mode:'lines+markers', name:'XFOIL', line:{color:'#1f77b4', width:2}, marker:{size:6}},
                {x:mf.alpha, y:ld(mf.cl,mf.cd), mode:'lines+markers', name:'Mfoil', line:{color:'#2ca02c', width:2, dash:'dash'}, marker:{size:6, symbol:'diamond'}},
                {x:rf.alpha, y:ld(rf.cl,rf.cd), mode:'lines+markers', name:'RustFoil', line:{color:'#ff7f0e', width:2, dash:'dot'}, marker:{size:6, symbol:'square'}}
            ], {...layout, title:'<b>L/D Ratio</b>', xaxis:{title:'α (deg)', gridcolor:'#eee'}, yaxis:{title:'L/D', gridcolor:'#eee', rangemode:'tozero'}});
        }
        updatePlots();
    </script>
</body>
</html>"""
    
    html = html.replace('DATA_PLACEHOLDER', json.dumps(all_data))
    
    with open(output_path, 'w') as f:
        f.write(html)
    
    print(f"\nSaved: {output_path}")


def main():
    print("=" * 60)
    print("Quick Three-Solver Comparison")
    print("=" * 60)
    
    # Load existing XFOIL/RustFoil data
    existing = load_existing_data()
    
    all_data = []
    
    for naca, reynolds in CONFIGS:
        name = f"NACA {naca}"
        re_key = f"Re{reynolds:.0e}"
        
        # Get existing data
        xfoil_data = None
        rustfoil_data = None
        
        if name in existing and re_key in existing[name]:
            xfoil_data = existing[name][re_key].get('xfoil')
            rustfoil_data = existing[name][re_key].get('rustfoil')
            print(f"\n{name} Re={reynolds:.0e}: loaded XFOIL/RustFoil from cache")
        
        # Run Mfoil
        mfoil_data = run_mfoil(naca, reynolds, ALPHAS)
        
        all_data.append({
            'airfoil': naca,
            'name': name,
            'reynolds': reynolds,
            'xfoil': xfoil_data,
            'mfoil': mfoil_data,
            'rustfoil': rustfoil_data,
        })
    
    # Generate HTML
    html_path = OUTPUT_DIR / "three_solver_comparison.html"
    generate_html(all_data, html_path)
    
    # Save JSON
    json_path = OUTPUT_DIR / "three_solver_data.json"
    with open(json_path, 'w') as f:
        json.dump(all_data, f, indent=2)
    print(f"Saved: {json_path}")
    
    print("\n" + "=" * 60)
    print(f"Open: file://{html_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()

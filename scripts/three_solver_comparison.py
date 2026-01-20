#!/usr/bin/env python3
"""
Three-Solver Comparison: XFOIL vs Mfoil vs RustFoil

Creates an interactive Plotly HTML plot comparing aerodynamic coefficients
from all three solvers for multiple airfoils and Reynolds numbers.
"""

import subprocess
import json
import os
import sys
from pathlib import Path
import tempfile
import re
import numpy as np

# Add mfoil to path
PROJECT_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_DIR / "mfoil"))

from mfoil import mfoil

# Configuration
DEFAULT_XFOIL = PROJECT_DIR / "Xfoil" / "bin" / "xfoil"
XFOIL_PATH = os.environ.get(
    "XFOIL_PATH",
    str(DEFAULT_XFOIL if DEFAULT_XFOIL.exists() else "/Applications/ESP128/EngSketchPad/bin/xfoil")
)
OUTPUT_DIR = PROJECT_DIR / "testdata" / "xfoil_comparison"

# Test matrix
AIRFOILS = [
    {"name": "NACA 0012", "naca": "0012"},
    {"name": "NACA 2412", "naca": "2412"},
    {"name": "NACA 4412", "naca": "4412"},
]

REYNOLDS = [1e6, 3e6]
NCRIT = 9.0
ALPHA_RANGE = list(range(-4, 15, 1))  # -4 to 14 degrees


def run_xfoil_polar(naca: str, reynolds: float, alphas: list, ncrit: float = 9.0) -> dict:
    """Run XFOIL and extract polar data."""
    
    if not Path(XFOIL_PATH).exists():
        print(f"  XFOIL not found at: {XFOIL_PATH}")
        return None

    alpha_start, alpha_end = min(alphas), max(alphas)
    alpha_step = 1

    with tempfile.TemporaryDirectory() as tmpdir:
        polar_filename = "polar.txt"
        polar_file = os.path.join(tmpdir, polar_filename)
        
        commands = f"""NACA {naca}
PANE
OPER
VISC {reynolds:.0f}
VPAR
N {ncrit}

PACC
{polar_filename}

ASEQ {alpha_start} {alpha_end} {alpha_step}

QUIT
"""
        
        try:
            result = subprocess.run(
                [XFOIL_PATH],
                input=commands,
                capture_output=True,
                text=True,
                timeout=120,
                cwd=tmpdir
            )
        except subprocess.TimeoutExpired:
            print(f"  XFOIL timeout for NACA {naca} at Re={reynolds:.0e}")
            return None
        except Exception as e:
            print(f"  XFOIL error for NACA {naca}: {e}")
            return None
        
        if not os.path.exists(polar_file):
            return None
        
        return parse_xfoil_polar(polar_file)


def parse_xfoil_polar(filepath: str) -> dict:
    """Parse XFOIL polar output file."""
    data = {"alpha": [], "cl": [], "cd": [], "cm": []}
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        if not line or not re.match(r'^\s*-?\d', line):
            continue
        
        parts = line.split()
        if len(parts) >= 5:
            try:
                alpha = float(parts[0])
                cl = float(parts[1])
                cd = float(parts[2])
                cm = float(parts[4])
                
                if abs(alpha) > 30 or abs(cl) > 5 or cd < 0 or cd > 1:
                    continue
                
                data["alpha"].append(alpha)
                data["cl"].append(cl)
                data["cd"].append(cd)
                data["cm"].append(cm)
            except (ValueError, IndexError):
                continue
    
    return data if data["alpha"] else None


def run_mfoil_polar(naca: str, reynolds: float, alphas: list) -> dict:
    """Run Mfoil and extract polar data."""
    import signal
    import warnings
    warnings.filterwarnings('ignore')
    
    data = {"alpha": [], "cl": [], "cd": [], "cm": []}
    
    try:
        m = mfoil(naca=naca, npanel=149)  # Fewer panels for speed
        m.param.verb = 0  # Quiet mode
        m.param.doplot = False
        m.param.niglob = 25  # Fewer iterations
        
        for alpha in alphas:
            try:
                m.setoper(alpha=alpha, Re=reynolds, Ma=0.0)
                m.oper.initbl = True  # Re-initialize BL for each alpha
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
            except Exception as e:
                # Skip this alpha on error
                pass
                
    except Exception as e:
        print(f"  Mfoil error for NACA {naca}: {e}")
        return None
    
    return data if data["alpha"] else None


def ensure_rustfoil_binary() -> Path:
    """Build or locate the RustFoil polar_gen binary."""
    solver_dir = PROJECT_DIR / "crates" / "rustfoil-solver"
    candidates = [
        PROJECT_DIR / "target" / "release" / "examples" / "polar_gen",
        solver_dir / "target" / "release" / "examples" / "polar_gen",
    ]
    for binary in candidates:
        if binary.exists():
            return binary

    print("  Building RustFoil polar_gen (release)...")
    result = subprocess.run(
        ["cargo", "build", "--release", "--example", "polar_gen"],
        capture_output=True,
        text=True,
        cwd=str(solver_dir),
    )
    if result.returncode != 0:
        print(f"  RustFoil build error: {result.stderr[:200]}")
        return None

    for binary in candidates:
        if binary.exists():
            return binary
    return None


def run_rustfoil_polar(naca: str, reynolds: float, alphas: list) -> dict:
    """Run RustFoil polar generation."""
    binary = ensure_rustfoil_binary()
    if not binary:
        return None

    alpha_start, alpha_end = min(alphas), max(alphas)
    alpha_step = 1

    try:
        result = subprocess.run(
            [str(binary), naca, str(reynolds), str(alpha_start), str(alpha_end), str(alpha_step)],
            capture_output=True,
            text=True,
            timeout=300,
            cwd=str(binary.parent),
        )
        
        if result.returncode != 0:
            print(f"  RustFoil error: {result.stderr[:200]}")
            return None
        
        output = result.stdout.strip()
        json_start = output.find('{')
        json_end = output.rfind('}') + 1
        if json_start >= 0 and json_end > json_start:
            json_str = output[json_start:json_end]
            return json.loads(json_str)
        
    except subprocess.TimeoutExpired:
        print(f"  RustFoil timeout for NACA {naca} at Re={reynolds:.0e}")
    except json.JSONDecodeError as e:
        print(f"  JSON parse error: {e}")
    except Exception as e:
        print(f"  RustFoil error: {e}")
    
    return None


def generate_html(all_data: list, output_path: Path):
    """Generate interactive HTML comparison plot."""
    
    html_template = """<!DOCTYPE html>
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
        .container {
            max-width: 1600px;
            margin: 0 auto;
        }
        h1 { 
            color: white; 
            text-align: center;
            margin-bottom: 10px;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
        }
        .subtitle {
            color: rgba(255,255,255,0.9);
            text-align: center;
            margin-bottom: 20px;
            font-size: 14px;
        }
        .info { 
            background: rgba(255,255,255,0.95); 
            padding: 15px 20px; 
            border-radius: 10px; 
            margin: 15px 0;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }
        .info h3 { margin: 0 0 10px 0; color: #333; }
        .info ul { margin: 0; padding-left: 20px; }
        .info li { margin: 3px 0; color: #555; }
        .selector { 
            background: rgba(255,255,255,0.95);
            padding: 15px 20px;
            border-radius: 10px;
            margin: 15px 0;
            display: flex;
            align-items: center;
            gap: 20px;
            flex-wrap: wrap;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }
        .selector label { font-weight: 600; color: #333; }
        select { 
            padding: 10px 15px; 
            font-size: 14px; 
            border: 2px solid #ddd;
            border-radius: 6px;
            background: white;
            cursor: pointer;
            min-width: 150px;
        }
        select:hover { border-color: #667eea; }
        .plot-container { 
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 20px;
            margin-top: 20px;
        }
        .plot { 
            background: white; 
            padding: 15px; 
            border-radius: 10px; 
            box-shadow: 0 4px 15px rgba(0,0,0,0.15);
        }
        .legend-info {
            background: rgba(255,255,255,0.95);
            padding: 10px 20px;
            border-radius: 10px;
            margin: 15px 0;
            display: flex;
            justify-content: center;
            gap: 30px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }
        .legend-item {
            display: flex;
            align-items: center;
            gap: 8px;
        }
        .legend-color {
            width: 30px;
            height: 4px;
            border-radius: 2px;
        }
        .xfoil-color { background: #1f77b4; }
        .mfoil-color { background: #2ca02c; }
        .rustfoil-color { background: #ff7f0e; }
        .stats-table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 10px;
            font-size: 13px;
        }
        .stats-table th, .stats-table td {
            padding: 8px 12px;
            text-align: left;
            border-bottom: 1px solid #eee;
        }
        .stats-table th { background: #f5f5f5; font-weight: 600; }
        @media (max-width: 1200px) {
            .plot-container { grid-template-columns: 1fr; }
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>Three-Solver Aerodynamic Comparison</h1>
        <p class="subtitle">XFOIL (Fortran) vs Mfoil (Python) vs RustFoil (Rust)</p>
        
        <div class="legend-info">
            <div class="legend-item">
                <div class="legend-color xfoil-color"></div>
                <span><strong>XFOIL</strong> - Reference (Drela)</span>
            </div>
            <div class="legend-item">
                <div class="legend-color mfoil-color"></div>
                <span><strong>Mfoil</strong> - Python (Fidkowski)</span>
            </div>
            <div class="legend-item">
                <div class="legend-color rustfoil-color"></div>
                <span><strong>RustFoil</strong> - This Project</span>
            </div>
        </div>
        
        <div class="selector">
            <label>Airfoil:</label>
            <select id="airfoilSelect" onchange="updatePlots()">
                <option value="0012">NACA 0012</option>
                <option value="2412">NACA 2412</option>
                <option value="4412">NACA 4412</option>
            </select>
            <label>Reynolds Number:</label>
            <select id="reSelect" onchange="updatePlots()">
                <option value="1000000">1 Million</option>
                <option value="3000000">3 Million</option>
            </select>
        </div>
        
        <div class="info" id="statsInfo">
            <h3>Comparison Statistics</h3>
            <table class="stats-table" id="statsTable">
                <tr>
                    <th>Metric</th>
                    <th>Mfoil vs XFOIL</th>
                    <th>RustFoil vs XFOIL</th>
                </tr>
            </table>
        </div>
        
        <div class="plot-container">
            <div class="plot" id="clPlot" style="height: 450px;"></div>
            <div class="plot" id="cdPlot" style="height: 450px;"></div>
            <div class="plot" id="polarPlot" style="height: 450px;"></div>
            <div class="plot" id="ldPlot" style="height: 450px;"></div>
        </div>
    </div>
    
    <script>
        const allData = ALL_DATA_PLACEHOLDER;
        
        function computeStats(ref, test) {
            if (!ref || !test || !ref.alpha || !test.alpha) return null;
            
            let clErrs = [], cdErrs = [];
            for (let i = 0; i < ref.alpha.length; i++) {
                const alpha = ref.alpha[i];
                const j = test.alpha.findIndex(a => Math.abs(a - alpha) < 0.5);
                if (j >= 0 && test.cl[j] !== null && test.cd[j] !== null) {
                    clErrs.push(test.cl[j] - ref.cl[i]);
                    cdErrs.push(test.cd[j] - ref.cd[i]);
                }
            }
            
            if (clErrs.length === 0) return null;
            
            const rms = arr => Math.sqrt(arr.reduce((s, x) => s + x*x, 0) / arr.length);
            const maxAbs = arr => Math.max(...arr.map(Math.abs));
            
            return {
                clRms: rms(clErrs),
                cdRms: rms(cdErrs),
                clMax: maxAbs(clErrs),
                cdMax: maxAbs(cdErrs),
                nPoints: clErrs.length
            };
        }
        
        function updateStats(xfoil, mfoil, rustfoil) {
            const mfoilStats = computeStats(xfoil, mfoil);
            const rustfoilStats = computeStats(xfoil, rustfoil);
            
            let html = `
                <tr><th>Metric</th><th>Mfoil vs XFOIL</th><th>RustFoil vs XFOIL</th></tr>
            `;
            
            if (mfoilStats && rustfoilStats) {
                html += `
                    <tr><td>Cl RMS Error</td><td>${mfoilStats.clRms.toFixed(4)}</td><td>${rustfoilStats.clRms.toFixed(4)}</td></tr>
                    <tr><td>Cd RMS Error (counts)</td><td>${(mfoilStats.cdRms * 10000).toFixed(1)}</td><td>${(rustfoilStats.cdRms * 10000).toFixed(1)}</td></tr>
                    <tr><td>Cl Max Error</td><td>${mfoilStats.clMax.toFixed(4)}</td><td>${rustfoilStats.clMax.toFixed(4)}</td></tr>
                    <tr><td>Cd Max Error (counts)</td><td>${(mfoilStats.cdMax * 10000).toFixed(1)}</td><td>${(rustfoilStats.cdMax * 10000).toFixed(1)}</td></tr>
                    <tr><td>Comparison Points</td><td>${mfoilStats.nPoints}</td><td>${rustfoilStats.nPoints}</td></tr>
                `;
            } else {
                html += `<tr><td colspan="3">Insufficient data for comparison</td></tr>`;
            }
            
            document.getElementById('statsTable').innerHTML = html;
        }
        
        function updatePlots() {
            const airfoil = document.getElementById('airfoilSelect').value;
            const re = parseFloat(document.getElementById('reSelect').value);
            
            const match = allData.find(d => d.airfoil === airfoil && d.reynolds === re);
            if (!match) {
                console.log('No data for', airfoil, re);
                return;
            }
            
            const xfoil = match.xfoil || {alpha: [], cl: [], cd: [], cm: []};
            const mfoil = match.mfoil || {alpha: [], cl: [], cd: [], cm: []};
            const rustfoil = match.rustfoil || {alpha: [], cl: [], cd: [], cm: []};
            
            updateStats(xfoil, mfoil, rustfoil);
            
            const commonLayout = {
                font: { family: '-apple-system, BlinkMacSystemFont, Segoe UI, Roboto, sans-serif' },
                legend: { x: 0.02, y: 0.98, bgcolor: 'rgba(255,255,255,0.8)' },
                margin: { t: 50, r: 30, b: 50, l: 60 }
            };
            
            // Cl vs Alpha
            Plotly.newPlot('clPlot', [
                { x: xfoil.alpha, y: xfoil.cl, mode: 'lines+markers', name: 'XFOIL', 
                  line: {color: '#1f77b4', width: 2}, marker: {size: 6} },
                { x: mfoil.alpha, y: mfoil.cl, mode: 'lines+markers', name: 'Mfoil',
                  line: {color: '#2ca02c', width: 2, dash: 'dash'}, marker: {size: 6, symbol: 'diamond'} },
                { x: rustfoil.alpha, y: rustfoil.cl, mode: 'lines+markers', name: 'RustFoil',
                  line: {color: '#ff7f0e', width: 2, dash: 'dot'}, marker: {size: 6, symbol: 'square'} }
            ], { 
                ...commonLayout,
                title: '<b>Lift Coefficient vs Angle of Attack</b>', 
                xaxis: {title: 'α (deg)', gridcolor: '#eee'}, 
                yaxis: {title: 'Cl', gridcolor: '#eee'}
            });
            
            // Cd vs Alpha
            Plotly.newPlot('cdPlot', [
                { x: xfoil.alpha, y: xfoil.cd.map(c => c * 10000), mode: 'lines+markers', name: 'XFOIL',
                  line: {color: '#1f77b4', width: 2}, marker: {size: 6} },
                { x: mfoil.alpha, y: mfoil.cd.map(c => c * 10000), mode: 'lines+markers', name: 'Mfoil',
                  line: {color: '#2ca02c', width: 2, dash: 'dash'}, marker: {size: 6, symbol: 'diamond'} },
                { x: rustfoil.alpha, y: rustfoil.cd.map(c => c * 10000), mode: 'lines+markers', name: 'RustFoil',
                  line: {color: '#ff7f0e', width: 2, dash: 'dot'}, marker: {size: 6, symbol: 'square'} }
            ], { 
                ...commonLayout,
                title: '<b>Drag Coefficient vs Angle of Attack</b>', 
                xaxis: {title: 'α (deg)', gridcolor: '#eee'}, 
                yaxis: {title: 'Cd (counts)', gridcolor: '#eee'}
            });
            
            // Drag Polar (Cl vs Cd)
            Plotly.newPlot('polarPlot', [
                { x: xfoil.cd.map(c => c * 10000), y: xfoil.cl, mode: 'lines+markers', name: 'XFOIL',
                  line: {color: '#1f77b4', width: 2}, marker: {size: 6} },
                { x: mfoil.cd.map(c => c * 10000), y: mfoil.cl, mode: 'lines+markers', name: 'Mfoil',
                  line: {color: '#2ca02c', width: 2, dash: 'dash'}, marker: {size: 6, symbol: 'diamond'} },
                { x: rustfoil.cd.map(c => c * 10000), y: rustfoil.cl, mode: 'lines+markers', name: 'RustFoil',
                  line: {color: '#ff7f0e', width: 2, dash: 'dot'}, marker: {size: 6, symbol: 'square'} }
            ], { 
                ...commonLayout,
                title: '<b>Drag Polar</b>', 
                xaxis: {title: 'Cd (counts)', gridcolor: '#eee'}, 
                yaxis: {title: 'Cl', gridcolor: '#eee'},
                legend: { x: 0.98, y: 0.02, xanchor: 'right', bgcolor: 'rgba(255,255,255,0.8)' }
            });
            
            // L/D vs Alpha
            const ld = (cl, cd) => cl.map((c, i) => cd[i] > 1e-6 ? c / cd[i] : null);
            Plotly.newPlot('ldPlot', [
                { x: xfoil.alpha, y: ld(xfoil.cl, xfoil.cd), mode: 'lines+markers', name: 'XFOIL',
                  line: {color: '#1f77b4', width: 2}, marker: {size: 6} },
                { x: mfoil.alpha, y: ld(mfoil.cl, mfoil.cd), mode: 'lines+markers', name: 'Mfoil',
                  line: {color: '#2ca02c', width: 2, dash: 'dash'}, marker: {size: 6, symbol: 'diamond'} },
                { x: rustfoil.alpha, y: ld(rustfoil.cl, rustfoil.cd), mode: 'lines+markers', name: 'RustFoil',
                  line: {color: '#ff7f0e', width: 2, dash: 'dot'}, marker: {size: 6, symbol: 'square'} }
            ], { 
                ...commonLayout,
                title: '<b>Lift-to-Drag Ratio</b>', 
                xaxis: {title: 'α (deg)', gridcolor: '#eee'}, 
                yaxis: {title: 'L/D', gridcolor: '#eee', rangemode: 'tozero'}
            });
        }
        
        // Initialize
        updatePlots();
    </script>
</body>
</html>
"""
    
    # Replace placeholder with actual data
    html_content = html_template.replace('ALL_DATA_PLACEHOLDER', json.dumps(all_data))
    
    with open(output_path, 'w') as f:
        f.write(html_content)
    
    print(f"\nSaved interactive comparison: {output_path}")


def main():
    print("=" * 70)
    print("Three-Solver Comparison: XFOIL vs Mfoil vs RustFoil")
    print("=" * 70)
    
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Check for cached data
    json_path = OUTPUT_DIR / "three_solver_data.json"
    cached_data = {}
    if json_path.exists():
        try:
            with open(json_path, 'r') as f:
                cached_list = json.load(f)
                for item in cached_list:
                    key = f"{item['airfoil']}_{item['reynolds']}"
                    cached_data[key] = item
            print(f"Loaded {len(cached_data)} cached results")
        except:
            pass
    
    all_data = []
    
    for airfoil in AIRFOILS:
        print(f"\n{airfoil['name']}")
        print("-" * 40)
        
        for reynolds in REYNOLDS:
            print(f"\n  Re = {reynolds:.0e}")
            cache_key = f"{airfoil['naca']}_{reynolds}"
            
            # Check cache
            if cache_key in cached_data:
                cached = cached_data[cache_key]
                if cached.get('xfoil') and cached.get('mfoil') and cached.get('rustfoil'):
                    print("    Using cached data ✓")
                    all_data.append(cached)
                    continue
            
            # Run XFOIL
            print("    Running XFOIL...", end=" ", flush=True)
            xfoil_data = run_xfoil_polar(airfoil['naca'], reynolds, ALPHA_RANGE, NCRIT)
            if xfoil_data:
                print(f"✓ ({len(xfoil_data['alpha'])} points)")
            else:
                print("✗")
            
            # Run Mfoil
            print("    Running Mfoil...", end=" ", flush=True)
            mfoil_data = run_mfoil_polar(airfoil['naca'], reynolds, ALPHA_RANGE)
            if mfoil_data:
                print(f"✓ ({len(mfoil_data['alpha'])} points)")
            else:
                print("✗")
            
            # Run RustFoil
            print("    Running RustFoil...", end=" ", flush=True)
            rustfoil_data = run_rustfoil_polar(airfoil['naca'], reynolds, ALPHA_RANGE)
            if rustfoil_data:
                print(f"✓ ({len(rustfoil_data['alpha'])} points)")
            else:
                print("✗")
            
            result = {
                'airfoil': airfoil['naca'],
                'name': airfoil['name'],
                'reynolds': reynolds,
                'xfoil': xfoil_data,
                'mfoil': mfoil_data,
                'rustfoil': rustfoil_data,
            }
            all_data.append(result)
            
            # Save incrementally
            with open(json_path, 'w') as f:
                json.dump(all_data, f, indent=2)
    
    # Generate interactive HTML
    html_path = OUTPUT_DIR / "three_solver_comparison.html"
    generate_html(all_data, html_path)
    
    # Also save JSON data
    json_path = OUTPUT_DIR / "three_solver_data.json"
    with open(json_path, 'w') as f:
        json.dump(all_data, f, indent=2)
    print(f"Saved JSON data: {json_path}")
    
    print("\n" + "=" * 70)
    print("Comparison complete!")
    print(f"Open in browser: file://{html_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()

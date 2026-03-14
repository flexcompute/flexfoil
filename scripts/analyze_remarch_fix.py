#!/usr/bin/env python3
"""
Analyze the re-march fix at α=15° to verify H growth at TE
"""
import json
from collections import Counter

trace_file = "crates/rustfoil-solver/traces/rustfoil_new/rustfoil_alpha_15.json"

with open(trace_file, "r") as f:
    data = json.load(f)

# Extract events
events = data.get('events', [])

# Count MRCHUE events per station (upper surface only)
mrchue_upper = [e for e in events if e.get('subroutine') == 'MRCHUE' and e.get('side') == 1]

# Group by ibl and get x values
ibl_x_map = {}
for e in mrchue_upper:
    ibl = e.get('ibl', 0)
    x = e.get('x', 0)
    if ibl not in ibl_x_map:
        ibl_x_map[ibl] = x

# Find the last station on the surface (x <= 1.0)
surface_ibls = [(ibl, x) for ibl, x in ibl_x_map.items() if x <= 1.01]
surface_ibls.sort(key=lambda t: t[1], reverse=True)

if surface_ibls:
    te_ibl, te_x = surface_ibls[0]
else:
    # Fallback: find station with x closest to 1.0
    te_ibl = max(ibl_x_map.keys(), key=lambda ibl: -abs(ibl_x_map[ibl] - 1.0))
    te_x = ibl_x_map[te_ibl]

# Count occurrences
ibl_counts = Counter(e.get('ibl', 0) for e in mrchue_upper)

print(f"=== Re-March Fix Verification at α=15° ===")
print(f"\n1. MRCHUE EVENT COUNT PER STATION (Upper Surface)")
print(f"   Total unique stations: {len(ibl_counts)}")
print(f"   Trailing edge station: ibl={te_ibl} (x≈{te_x:.4f})")
print(f"   TE appears {ibl_counts[te_ibl]} times (should be ~8-20 if re-march works)")

# Track H evolution at TE - look for the station closest to x=1.0
te_events = [e for e in mrchue_upper if e.get('ibl') == te_ibl]

print(f"\n2. H EVOLUTION AT TE (ibl={te_ibl}, x≈{te_x:.4f}):")
if te_events:
    print(f"   {'Occurrence':<12} {'H (Hk)':<12} {'Cf':<12} {'Ue':<12} {'theta':<12} {'dstar':<12}")
    print(f"   {'-'*78}")
    
    for i, e in enumerate(te_events):
        h = e.get('Hk', 0)
        cf = e.get('Cf', 0)
        ue = e.get('Ue', 0)
        theta = e.get('theta', 0)
        dstar = e.get('delta_star', 0)
        print(f"   {i+1:<12} {h:<12.3f} {cf:<12.6f} {ue:<12.4f} {theta:<12.6f} {dstar:<12.6f}")
    
    # Find max H at TE
    max_h = max(e.get('Hk', 0) for e in te_events)
    min_h = min(e.get('Hk', 0) for e in te_events)
    first_h = te_events[0].get('Hk', 0)
    last_h = te_events[-1].get('Hk', 0)
    
    print(f"\n3. H STATISTICS AT TE:")
    print(f"   First H:  {first_h:.3f}")
    print(f"   Last H:   {last_h:.3f}")
    print(f"   Min H:    {min_h:.3f}")
    print(f"   Max H:    {max_h:.3f}")
    print(f"   Range:    {max_h - min_h:.3f}")
    print(f"   Growth:   {last_h - first_h:.3f} ({100*(last_h - first_h)/first_h:.1f}%)")
    
    # Check success criteria
    print(f"\n4. SUCCESS CRITERIA:")
    print(f"   ✓ PASS if ibl={te_ibl} appears 8-20 times: {ibl_counts[te_ibl]} times - {'✅ PASS' if 8 <= ibl_counts[te_ibl] <= 20 else '❌ FAIL'}")
    print(f"   ✓ PASS if H grows from ~2.4 → ~3.0: {first_h:.3f} → {last_h:.3f}")
    
    if first_h < 2.5 and last_h > 2.8:
        print(f"      {'✅ PASS'} - H showed strong growth")
    elif last_h > 2.5:
        print(f"      {'⚠️  MARGINAL'} - H grew but didn't reach 3.0")
    else:
        print(f"      {'❌ FAIL'} - H didn't grow sufficiently")
    
    print(f"   ✓ PASS if Max H ≥ 2.9: {max_h:.3f} - {'✅ PASS' if max_h >= 2.9 else '⚠️  MARGINAL' if max_h >= 2.5 else '❌ FAIL'}")
else:
    print("   No events found!")

# Show station coverage near TE
print(f"\n5. STATION COVERAGE (showing last 10 stations):")
sorted_stations = sorted(ibl_counts.items(), reverse=True)[:10]
for ibl, count in sorted_stations:
    # Find x position
    x_pos = ibl_x_map.get(ibl, 0)
    print(f"   ibl={ibl:<3} (x≈{x_pos:.4f}): {count} occurrences")

# Check if we actually have multiple marches per iteration
print(f"\n6. RE-MARCH VERIFICATION:")
print(f"   If re-march is working, each station should appear once per global iteration.")
print(f"   Station ibl={te_ibl} appears {ibl_counts[te_ibl]} times")

# Let's look at a few other surface stations
print(f"\n7. OTHER SURFACE STATIONS (for comparison):")
sample_stations = [te_ibl - 10, te_ibl - 5, te_ibl - 1] if te_ibl > 10 else []
for ibl in sample_stations:
    if ibl in ibl_counts and ibl in ibl_x_map:
        count = ibl_counts[ibl]
        x = ibl_x_map[ibl]
        print(f"   ibl={ibl:<3} (x≈{x:.4f}): {count} occurrences")

#!/usr/bin/env python3
"""Compare RustFoil and XFOIL polar results."""

# XFOIL results (from full polar)
xfoil_data = """
-4, -0.4278, 0.00728, 0.9685, 0.2537
-3, -0.3199, 0.00639, 0.9285, 0.3642
-2, -0.2142, 0.00580, 0.8676, 0.4743
-1, -0.1074, 0.00549, 0.7848, 0.5826
0, 0.0000, 0.00541, 0.6870, 0.6870
1, 0.1074, 0.00549, 0.5826, 0.7848
2, 0.2142, 0.00580, 0.4743, 0.8676
3, 0.3200, 0.00639, 0.3642, 0.9285
4, 0.4278, 0.00728, 0.2537, 0.9685
5, 0.5580, 0.00848, 0.1485, 0.9849
6, 0.6948, 0.00973, 0.0812, 0.9940
7, 0.8264, 0.01094, 0.0506, 1.0000
8, 0.9098, 0.01211, 0.0381, 1.0000
9, 0.9947, 0.01341, 0.0307, 1.0000
10, 1.0809, 0.01498, 0.0256, 1.0000
"""

# RustFoil results (from viscous-polar output)
rustfoil_data = """
-4, -0.4677, 0.01728, 0.9619, 0.0062
-3, -0.3508, 0.01296, 0.9263, 0.0062
-2, -0.2339, 0.01261, 0.8802, 0.0062
-1, -0.1170, 0.00964, 0.8247, 0.6167
0, 0.0000, 0.01672, 0.0381, 0.7612
1, 0.1170, 0.00964, 0.6167, 0.8247
2, 0.2339, 0.01261, 0.0062, 0.8802
3, 0.3508, 0.01296, 0.0062, 0.9263
4, 0.4677, 0.01728, 0.0062, 0.9619
5, 0.5844, 0.04731, 0.0015, 0.9862
6, 0.7011, 0.10274, 0.0015, 1.0000
7, 0.8176, 0.12194, 0.0015, 1.0000
8, 0.9339, 0.16012, 0.0015, 1.0000
9, 1.0501, 0.13924, 0.0015, 1.0000
10, 1.1661, 0.10929, 0.0015, 1.0000
"""

def parse_data(data_str):
    result = []
    for line in data_str.strip().split('\n'):
        if line.strip():
            parts = [float(x.strip()) for x in line.split(',')]
            result.append(parts)
    return result

xfoil = parse_data(xfoil_data)
rustfoil = parse_data(rustfoil_data)

print("=" * 90)
print("POLAR COMPARISON: RustFoil vs XFOIL (NACA 0012, Re=1e6)")
print("=" * 90)
print()
print(f"{'Alpha':>6} | {'CL (RF)':>8} {'CL (XF)':>8} {'Diff':>7} | {'CD (RF)':>8} {'CD (XF)':>8} {'Ratio':>7} | {'xtr_u RF':>8} {'xtr_u XF':>8}")
print("-" * 90)

for rf, xf in zip(rustfoil, xfoil):
    alpha = rf[0]
    cl_rf, cd_rf, xtr_u_rf = rf[1], rf[2], rf[3]
    cl_xf, cd_xf, xtr_u_xf = xf[1], xf[2], xf[3]
    
    cl_diff = cl_rf - cl_xf
    cd_ratio = cd_rf / cd_xf if cd_xf > 0 else 0
    
    print(f"{alpha:6.1f} | {cl_rf:8.4f} {cl_xf:8.4f} {cl_diff:+7.4f} | {cd_rf:8.5f} {cd_xf:8.5f} {cd_ratio:7.2f}x | {xtr_u_rf:8.4f} {xtr_u_xf:8.4f}")

print("-" * 90)
print()
print("Key observations:")
print("- CL values: RustFoil inviscid CL (used directly), XFOIL viscous CL (slightly lower)")
print("- CD values: RustFoil CD is 1.5-10x higher than XFOIL, especially at α > 5°")
print("- x_tr: RustFoil showing very early transition (0.0015-0.0062) at higher α - suspicious!")

declare module 'rustfoil-wasm' {
  const init: () => Promise<void>;
  export default init;

  export const generate_naca4: (...args: any[]) => any;
  export const generate_naca4_from_string: (...args: any[]) => any;
  export const generate_naca4_xfoil: (...args: any[]) => any;
  export const repanel_with_spacing_and_curvature: (...args: any[]) => any;
  export const repanel_xfoil: (...args: any[]) => any;
  export const repanel_xfoil_with_params: (...args: any[]) => any;
  export const compute_curvature_spacing: (...args: any[]) => any;
  export const analyze_airfoil: (...args: any[]) => any;
  export const analyze_airfoil_faithful: (...args: any[]) => any;
  export const compute_streamlines: (...args: any[]) => any;
  export const compute_streamlines_faithful: (...args: any[]) => any;
  export const compute_dividing_streamline: (...args: any[]) => any;
  export const compute_dividing_streamline_faithful: (...args: any[]) => any;
  export const compute_psi_grid: (...args: any[]) => any;
  export const compute_psi_grid_faithful: (...args: any[]) => any;
  export const get_bl_distribution_faithful: (...args: any[]) => any;
  export const get_bl_visualization_faithful: (...args: any[]) => any;
  export const greet: (...args: any[]) => any;
  export const RustFoil: any;
  export class WasmSmokeSystem {
    constructor(...args: any[]);
    [key: string]: any;
  }
  export const lerp_points: (...args: any[]) => any;
  export const lerp_array: (...args: any[]) => any;
  export const lerp_streamlines: (...args: any[]) => any;
  export const lerp_morph_state: (...args: any[]) => any;
}
